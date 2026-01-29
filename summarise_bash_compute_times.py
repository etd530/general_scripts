#!/usr/bin/env python3


#### LIBS ####
from datetime import timedelta

times = ["0:43.22", "0:44.15", "0:36.05", "0:52.82", "1:00.31", "0:43.47", "0:55.66", "0:42.34", "0:39.83", "0:40.79", "0:52.19", "1:10.06", "0:47.20", "0:32.20", "0:37.69", "0:51.32", "0:47.70", "0:37.50", "0:33.36", "0:46.76", "1:06.49", "0:23.38", "0:35.79", "0:57.01", "0:58.36", "0:41.82", "0:34.78", "0:24.46", "0:29.18", "0:46.20", "0:23.93", "0:35.29", "0:28.74", "0:17.62", "0:35.10", "0:33.70", "0:40.09", "0:38.24", "0:53.64", "0:29.94", "0:28.71", "0:30.50", "0:36.49", "0:17.96", "0:24.54", "0:31.10", "0:33.20", "0:43.36", "0:26.69", "0:28.54", "0:28.22"]

#### FUNS ####
def parse_time(t):
    parts = t.split(":")
    if len(parts) == 3:
        h, m, s = parts
        return int(h)*3600 + int(m)*60 + float(s)
    else:
        m, s = parts
        return int(m)*60 + float(s)

def format_seconds(sec):
    h = int(sec // 3600)
    m = int((sec % 3600) // 60)
    s = sec % 60
    if h > 0:
        return f"{h}:{m:02d}:{s:05.2f}"
    else:
        return f"{m}:{s:05.2f}"

#### MAIN ####
seconds = [parse_time(t) for t in times]
avg = sum(seconds)/len(seconds)
max = max(seconds)
print(f"Average: {round(avg)} seconds ({format_seconds(avg)})")
print(f"Max: {round(max)} seconds ({format_seconds(max)})")
