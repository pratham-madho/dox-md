import os, re, sys, datetime

out = sys.argv[1] if len(sys.argv) > 1 else "prod_B.out"
total_ns = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
size = os.path.getsize(out)
with open(out, "rb") as f:
    f.seek(max(0, size - 20000))                       # the progress line is rewritten at the end of the file
    txt = f.read().decode(errors="ignore").replace("\r", "\n")
m = re.findall(r"step (\d+), will finish (.+?\d{4})", txt)
if not m:
    sys.exit("no progress line yet - give it a minute and try again")
step, fin = m[-1]
ns = int(step) * 2e-6                                  # 2 fs per step
end = datetime.datetime.strptime(" ".join(fin.split()), "%a %b %d %H:%M:%S %Y")
h = (end - datetime.datetime.now()).total_seconds() / 3600
print(f"{ns:.1f} ns of {total_ns:.0f} done ({100*ns/total_ns:.0f}%) | finish {end:%a %d %b %H:%M} | "
      f"{int(h)} h {int((h % 1) * 60)} min left (about {h/24:.1f} days)")
