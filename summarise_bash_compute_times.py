#!/usr/bin/env python3


#### LIBS ####
from datetime import timedelta

times = ["0:14.55", "0:14.51", "0:14.74", "0:14.67", "0:14.47", "0:14.73", "0:02.90", "0:03.28", "0:02.86", "0:02.98", "0:02.76", "0:02.49"]

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
