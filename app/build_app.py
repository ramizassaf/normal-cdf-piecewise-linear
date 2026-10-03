"""Inject optimised model constants from results/results.json into the app.
Run after compute_results.py:  python app/build_app.py  ->  docs/index.html"""
import json, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
r = json.load(open(os.path.join(ROOT, "results", "results.json")))
models = {k: {"knots": r[f"model_{k}"]["knots"], "values": r[f"model_{k}"]["values"]} for k in "DEF"}
tpl = open(os.path.join(ROOT, "app", "template.html"), encoding="utf-8").read()
out = tpl.replace("__MODELS__", json.dumps(models))
os.makedirs(os.path.join(ROOT, "docs"), exist_ok=True)
open(os.path.join(ROOT, "docs", "index.html"), "w", encoding="utf-8").write(out)
print("wrote docs/index.html")
