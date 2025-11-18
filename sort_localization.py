import json

if __name__ == "__main__":
    with open("pibs/ui/localization.json", "r", encoding="utf-8") as f:
        v = json.load(f)

    with open("pibs/ui/localization.json", "w", encoding="utf-8") as f:
        json.dump(v, f, sort_keys=True, ensure_ascii=False, indent="\t")
