from ruamel.yaml import YAML

yaml = YAML()

import argparse
import datetime

parser = argparse.ArgumentParser()

parser.add_argument("citation_cff_file")
parser.add_argument("--version")
parser.add_argument("--doi")
parser.add_argument(
    "--date_released",
    type=datetime.date.fromisoformat,
    default=datetime.date.today(),
    help="Date in iso format (e.g. 2025-12-26)",
)

args = parser.parse_args()
with open(args.citation_cff_file) as f:
    d = yaml.load(f)

d["version"] = args.version
d["identifiers"][0]["value"] = args.doi
d["date-released"] = args.date_released.isoformat()

with open(args.citation_cff_file, "w") as f:
    yaml.dump(d, f)
