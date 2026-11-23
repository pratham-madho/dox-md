function out(name, arr, n,   i) {
  printf "[ %s ]\n", name
  for (i = 1; i <= n; i++) { printf "%6d", arr[i]; if (i % 15 == 0) printf "\n" }
  if (n % 15 != 0) printf "\n"
}
BEGIN { memb_re = "^(" memb ")$" }
NR == 2 { N = $1 }
NR > 2 && NR <= N + 2 {
  i = NR - 2
  rn = substr($0, 6, 5); an = substr($0, 11, 5)
  gsub(/ /, "", rn); gsub(/ /, "", an)
  all[++na] = i
  if (rn ~ memb_re) { m[++nm] = i; if (rn == "DPPC" && an == "P") p[++np] = i } else { s[++ns] = i }
  if (rn == drug) d[++nd] = i
  if (rn == ion)  c[++nc] = i
}
END {
  out("System", all, na); out("MEMB", m, nm); out("SOLV", s, ns)
  out("DOX", d, nd);      out("CLA", c, nc);  out("DPPC_P", p, np)
}
