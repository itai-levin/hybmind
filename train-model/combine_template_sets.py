import pandas as pd
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Return a combined json with all the input template sets"
    )
    parser.add_argument("--template-set-paths", nargs="+", dest="template_set_paths")
    parser.add_argument("--output-prefix", dest="out_prefix", default="data/processed/")
    parser.add_argument("--output-name", dest="out_name", default="retro.templates.combined.json.gz")
    return parser.parse_args()


def combine_template_sets(templates):
    df = pd.read_json(templates[0])
    for i in templates[1:]:
        df = df.append(pd.read_json(i), ignore_index=True)
    return df


if __name__ == "__main__":
    args = parse_arguments()
    print(args.template_set_paths)
    df = combine_template_sets(args.template_set_paths)

    df.to_json(args.out_prefix + args.out_name, orient="records")
