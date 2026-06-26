#!/usr/bin/env python3
"""
ray_terrain_intersection.py

Wild boar ground-coordinate extraction via ray-DEM intersection.

A single consolidated ethogram drives everything through its `type` column:

  type = point     (Feb-Mar dates)  -> one coordinate per observation, taken at
                                        `time`; camera azimuth from the absolute
                                        gimbal yaw `gimbal_heading(degrees)`.
  type = interval  (April dates)     -> traveling sampled every
                                        SAMPLE_TRAVELING_SEC s from `time_start`
                                        (exclusive of `time_end`); all other
                                        behaviors yield one coordinate at
                                        `time_start`; camera azimuth from the
                                        aircraft `compass_heading(degrees)`.

Both azimuth sources were validated against UAV video for their respective
sessions (see README, "Camera heading").

INPUT
  data/Ethogram.csv             consolidated ethogram (.csv or .xlsx)
  data/flight_logs/<date>/*-Flight-Airdata.csv
  data/DEM_study_site.tif                       1-m LiDAR DEM (EPSG:32652)

OUTPUT
  outputs/boar_coordinates.csv
  columns: date, time_kst, file, behavior, trackID, type, boar_lat, boar_lon
           (boar_lat/lon in WGS84, EPSG:4326)

ALTITUDE / CRS
  UAV altitude: flight-log MSL field `altitude_above_seaLevel(feet)` (N_M = 0.0).
  DEM CRS is read from the raster (dem.crs); DEM provided here is EPSG:32652.
  Output coordinates are WGS84 (EPSG:4326). See README/Methods for the
  vertical-reference rationale.

PATHS
  Resolved relative to the repository root by default. Each path can be
  overridden with an environment variable: WB_BASE_DIR, WB_ETHOGRAM,
  WB_FLIGHT_DIR, WB_DEM, WB_OUT.
"""

import os
import glob
from math import radians, sin, cos
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer

# ============================== CONFIG ==============================
BASE_DIR      = Path(os.environ.get("WB_BASE_DIR", "."))
DATA_DIR      = BASE_DIR / "data"
ETHOGRAM_PATH = Path(os.environ.get("WB_ETHOGRAM", DATA_DIR / "Ethogram.csv"))
FLIGHT_DIR    = Path(os.environ.get("WB_FLIGHT_DIR", DATA_DIR / "flight_logs"))
DEM_PATH      = Path(os.environ.get("WB_DEM", DATA_DIR / "DEM_study_site.tif"))
OUT_PATH      = Path(os.environ.get("WB_OUT", BASE_DIR / "outputs" / "boar_coordinates.csv"))

FLIGHT_CSV_PATTERN   = "*-Flight-Airdata.csv"
YEAR                 = 2025

HEADING_BY_TYPE = {                       # camera azimuth source per schema
    "point":    "gimbal_heading(degrees)",
    "interval": "compass_heading(degrees)",
}

N_M                  = 0.0      # geoid correction (m); 0 = shared vertical reference
SAMPLE_TRAVELING_SEC = 30       # interval sampling, traveling only
MATCH_TOL_SEC        = 10       # max |obs - flight-log| time gap for a match
MIN_GIMBAL_PITCH_DEG = -5.0     # drop near-horizontal / upward gimbal rows
RAY_STEP_COARSE      = 8.0
RAY_STEP_FINE        = 1.0
MAX_RANGE_M          = 3000.0
YAW_OFFSET_DEG       = 0.0
PITCH_OFFSET_DEG     = 0.0
TRAVELING_LABEL      = "traveling"
# ===================================================================


def _parse_hms(hms):
    if pd.isna(hms):
        return None
    s = str(hms).strip().replace("\u200b", "")
    for fmt in ("%H:%M:%S", "%H:%M:%S.%f"):
        t = pd.to_datetime(s, format=fmt, errors="coerce")
        if pd.notna(t):
            return t.hour * 3600 + t.minute * 60 + t.second + t.microsecond / 1e6
    t = pd.to_datetime(s, errors="coerce")
    if pd.notna(t):
        return t.hour * 3600 + t.minute * 60 + t.second + t.microsecond / 1e6
    # fallback: an Excel time cell stored as a day fraction (e.g. 0.857986 = 20:35:30)
    try:
        f = float(s)
        if 0 <= f < 1:
            return f * 86400.0
    except ValueError:
        pass
    return None


def _to_kst(base, hms):
    sec = _parse_hms(hms)
    return (base + pd.to_timedelta(sec, "s")) if sec is not None else pd.NaT


def fmt_date(date_str):
    return f"20{date_str[:2]}-{date_str[2:4]}-{date_str[4:6]}"


def fmt_time(t):
    return "" if pd.isna(t) else f"{t.hour}:{t.strftime('%M:%S')}"


def load_ethogram(path):
    path = Path(path)
    eth = pd.read_excel(path) if path.suffix.lower() in (".xlsx", ".xls") else pd.read_csv(path)
    eth.columns = eth.columns.str.strip()
    for col in ("Behavior", "type"):
        if col in eth.columns:
            eth[col] = eth[col].astype(str).str.strip()
    eth["date"] = eth.apply(
        lambda r: f"{str(YEAR)[2:]}{int(r['Month']):02d}{int(r['Day']):02d}", axis=1)
    return eth


def read_flight_csv(path, heading_col, n_m=0.0):
    df = pd.read_csv(path, low_memory=False)
    df.columns = df.columns.str.strip()
    ts_kst = (pd.to_datetime(df["datetime(utc)"], errors="coerce", utc=True)
                .dt.tz_convert("Asia/Seoul").dt.tz_localize(None))
    keep = pd.DataFrame({
        "time_kst":       ts_kst,
        "lat":            pd.to_numeric(df["latitude"],  errors="coerce"),
        "lon":            pd.to_numeric(df["longitude"], errors="coerce"),
        "abs_alt_m":      pd.to_numeric(df["altitude_above_seaLevel(feet)"],
                                        errors="coerce") * 0.3048 + n_m,
        "gimbal_heading": pd.to_numeric(df[heading_col], errors="coerce"),
        "gimbal_pitch":   pd.to_numeric(df["gimbal_pitch(degrees)"], errors="coerce"),
    })
    keep = keep.dropna(subset=["time_kst", "lat", "lon", "abs_alt_m",
                               "gimbal_heading", "gimbal_pitch"])
    keep = keep[keep["gimbal_pitch"] < MIN_GIMBAL_PITCH_DEG]
    return keep.sort_values("time_kst").reset_index(drop=True)


def load_all_flights(paths, heading_col, n_m=0.0):
    out = []
    for p in paths:
        try:
            out.append(read_flight_csv(p, heading_col, n_m=n_m))
        except Exception as e:
            print(f"  [WARNING] {os.path.basename(p)}: {e}")
    if not out:
        return pd.DataFrame(columns=["time_kst", "lat", "lon", "abs_alt_m",
                                     "gimbal_heading", "gimbal_pitch"])
    return pd.concat(out, ignore_index=True).sort_values("time_kst").reset_index(drop=True)


def build_requests(g, type_val):
    """Expand one date's observations into per-coordinate requests."""
    reqs = []
    for _, r in g.iterrows():
        base = pd.Timestamp(year=YEAR, month=int(r["Month"]), day=int(r["Day"]))
        beh  = str(r["Behavior"]).strip()
        common = {"file": str(r.get("File", "")).strip(), "behavior": beh,
                  "trackID": str(r.get("TrackID", "")).strip(), "type": type_val}
        if type_val == "interval":
            t0 = _to_kst(base, r.get("time_start"))
            t1 = _to_kst(base, r.get("time_end"))
            if pd.isna(t0):
                continue
            if beh.lower() == TRAVELING_LABEL.lower() and pd.notna(t1) and t1 > t0:
                t = t0
                while t < t1:
                    reqs.append({"t": t, **common})
                    t += pd.Timedelta(seconds=SAMPLE_TRAVELING_SEC)
            else:
                reqs.append({"t": t0, **common})
        else:  # point
            t = _to_kst(base, r.get("time"))
            if pd.isna(t):
                continue
            reqs.append({"t": t, **common})
    return reqs


def ray_dem_intersect(lat0, lon0, z0_m, heading_deg, pitch_deg, dem, to_dem, band,
                      coarse_step=8.0, fine_step=1.0, max_range_m=3000.0,
                      yaw_offset=0.0, pitch_offset=0.0):
    if any(pd.isna(v) for v in [lat0, lon0, z0_m, heading_deg, pitch_deg]):
        return None
    pit = radians(float(pitch_deg) + float(pitch_offset))
    if pit >= 0:
        return None
    hdg = radians(float(heading_deg) + float(yaw_offset))
    vx_e = sin(hdg) * cos(pit); vy_n = cos(hdg) * cos(pit); vz = sin(pit)
    H, W = dem.height, dem.width
    last_h = None; s = 0.0
    while s <= max_range_m:
        lat = lat0 + vy_n * s / 111000.0
        lon = lon0 + vx_e * s / (111000.0 * cos(radians(lat0)))
        z   = float(z0_m) + vz * s
        xd, yd = to_dem.transform(lon, lat)
        r, c   = dem.index(xd, yd)
        if not (0 <= r < H and 0 <= c < W):
            break
        zdv = band[r, c]
        if np.ma.is_masked(zdv):
            s += coarse_step; continue
        h = z - float(zdv)
        if last_h is not None and h <= 0 < last_h:
            sf0 = max(0.0, s - coarse_step * 2); sf1 = min(max_range_m, s + coarse_step)
            lh2 = None; ll2 = None; sf = sf0
            while sf <= sf1:
                lat2 = lat0 + vy_n * sf / 111000.0
                lon2 = lon0 + vx_e * sf / (111000.0 * cos(radians(lat0)))
                z2   = float(z0_m) + vz * sf
                x2, y2 = to_dem.transform(lon2, lat2)
                r2, c2 = dem.index(x2, y2)
                if not (0 <= r2 < H and 0 <= c2 < W):
                    break
                zdv2 = band[r2, c2]
                if np.ma.is_masked(zdv2):
                    sf += fine_step; continue
                h2 = z2 - float(zdv2)
                if lh2 is not None and h2 <= 0 < lh2:
                    alpha = lh2 / (lh2 - h2 + 1e-9)
                    return (ll2[0] + alpha * (lat2 - ll2[0]),
                            ll2[1] + alpha * (lon2 - ll2[1]))
                lh2 = h2; ll2 = (lat2, lon2); sf += fine_step
            break
        last_h = h; s += coarse_step
    return None


def main():
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    eth = load_ethogram(ETHOGRAM_PATH)
    all_rows = []

    with rasterio.open(DEM_PATH) as dem:
        to_dem = Transformer.from_crs("EPSG:4326", dem.crs, always_xy=True)
        band   = dem.read(1, masked=True)

        for date_str, g in eth.groupby("date"):
            type_val    = str(g["type"].iloc[0]).strip()
            heading_col = HEADING_BY_TYPE.get(type_val)
            if heading_col is None:
                print(f"[WARNING] {date_str}: unknown type '{type_val}', skipped"); continue

            date_fmt = fmt_date(date_str)
            reqs     = build_requests(g, type_val)
            csvs     = sorted(glob.glob(str(FLIGHT_DIR / date_str / FLIGHT_CSV_PATTERN)))
            flight   = load_all_flights(csvs, heading_col, n_m=N_M)
            n_t = sum(1 for q in reqs if q["behavior"].lower() == TRAVELING_LABEL.lower())
            print(f"\n{date_str} ({date_fmt})  type={type_val}  heading={heading_col}"
                  f"\n  requests={len(reqs)} (traveling {n_t})  flight rows={len(flight)}")

            miss_match = miss_ray = n_valid = 0
            for q in reqs:
                t = q["t"]; boar_lat = boar_lon = np.nan
                if len(flight):
                    i = flight["time_kst"].searchsorted(t)
                    cands = [flight.iloc[j] for j in [i - 1, i] if 0 <= j < len(flight)]
                    flt = None
                    if cands:
                        best = min(cands, key=lambda x: abs((x["time_kst"] - t).total_seconds()))
                        if abs((best["time_kst"] - t).total_seconds()) <= MATCH_TOL_SEC:
                            flt = best
                    if flt is None:
                        miss_match += 1
                    else:
                        hit = ray_dem_intersect(
                            float(flt["lat"]), float(flt["lon"]), float(flt["abs_alt_m"]),
                            float(flt["gimbal_heading"]), float(flt["gimbal_pitch"]),
                            dem, to_dem, band,
                            coarse_step=RAY_STEP_COARSE, fine_step=RAY_STEP_FINE,
                            max_range_m=MAX_RANGE_M,
                            yaw_offset=YAW_OFFSET_DEG, pitch_offset=PITCH_OFFSET_DEG)
                        if hit:
                            boar_lat, boar_lon = hit; n_valid += 1
                        else:
                            miss_ray += 1
                all_rows.append({
                    "date": date_fmt, "time_kst": fmt_time(t),
                    "file": q["file"], "behavior": q["behavior"],
                    "trackID": q["trackID"], "type": q["type"],
                    "boar_lat": boar_lat, "boar_lon": boar_lon,
                })
            print(f"  coords={n_valid}  no match={miss_match}  no ray={miss_ray}")

    out = pd.DataFrame(all_rows, columns=[
        "date", "time_kst", "file", "behavior", "trackID", "type", "boar_lat", "boar_lon"])
    out.to_csv(OUT_PATH, index=False)
    print(f"\nSaved: {OUT_PATH}  | rows={len(out)}  with coords={out['boar_lat'].notna().sum()}")


if __name__ == "__main__":
    main()
