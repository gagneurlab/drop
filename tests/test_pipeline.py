from tests.conftest import *
cores = "2"


def test_dryrun(demo_dir):
    r = run(["snakemake", "-n"], dir_path=demo_dir)
    message = "This was a dry-run (flag -n). The order of jobs does not reflect the order of execution."
    assert message in r.stdout


def test_dependencyGraph(demo_dir):
    r = run(["snakemake", "dependencyGraph", "-j"], dir_path=demo_dir)
    assert "aberrantExpression_dependency" in r.stderr
    assert "aberrantSplicing_dependency" in r.stderr
    assert "mae_dependency" in r.stderr
    assert "Finished job 0." in r.stderr
    # clean HTML output
    r = run(["snakemake", "clean", "-j"], dir_path=demo_dir)
    assert "clean" in r.stderr
    assert "Finished job 0." in r.stderr


def test_aberrantExpression(demo_dir):
    r = run(["snakemake", "aberrantExpression", "-j", cores], demo_dir)
    assert "Finished job 0." in r.stderr

    # check counts
    cnt_file = "Output/processed_data/aberrant_expression/v29/outrider/outrider/total_counts.Rds"
    r_cmd = """
        counts <- readRDS("{}")
        print(counts)
        """.format(cnt_file)
    r = runR(r_cmd, demo_dir)
    assert "dim: 805 12" in r.stdout
    assert "rownames(805): ENSG00000141956.13_3 ENSG00000141959.16_1 ..." in r.stdout

    # check OUTRIDER output
    output_dir = "Output/processed_results/aberrant_expression/v29/outrider/outrider"
    r_cmd = """ 
            # ods object
            ods <- readRDS(file.path("{}", "ods.Rds"))
            print(ods)
            
            # results table
            res <- readRDS(file.path("{}", "OUTRIDER_results_all.Rds"))
            print(paste(c("res:", dim(res)), collapse=" "))
            """.format(output_dir, output_dir)
    r = runR(r_cmd, demo_dir)
    assert "class: OutriderDataSet" in r.stdout
    assert "dim: 441 12" in r.stdout
    assert "res: 5292 15" in r.stdout


def test_aberrantSplicing(demo_dir):
    r = run(["snakemake", "aberrantSplicing", "-j", cores], demo_dir)
    assert "Finished job 0." in r.stderr

    # counts
    cnt_file = "Output/processed_data/aberrant_splicing/datasets/savedObjects/raw-fraser/fds-object.RDS"
    r_cmd = """ 
        library(FRASER)
        fds <- loadFraserDataSet(file="{}")
        print(fds)
        """.format(cnt_file)
    r = runR(r_cmd, demo_dir)
    assert "Number of samples:      10" in r.stdout
    assert "Number of junctions:    81" in r.stdout
    assert "Number of splice sites: 9" in r.stdout

    # FRASER results
    results_dir = "Output/processed_data/aberrant_splicing/results"
    r = run(f"wc -l {results_dir}/fraser_results_per_junction.tsv", demo_dir)
    assert "88 " in r.stdout
    r = run(f"wc -l {results_dir}/fraser_results.tsv", demo_dir)
    assert "1 " in r.stdout


def test_MAE(demo_dir):
    r = run(["snakemake", "mae", "-j", cores], demo_dir)
    assert "Finished job 0." in r.stderr

    # counts
    cnt_file = "Output/processed_data/mae/allelic_counts/HG00103--HG00103.4.M_120208_3.csv.gz"
    r_cmd = """ 
            library(data.table)
            cnts <- fread("{}")
            print(nrow(cnts))
            """.format(cnt_file)
    r = runR(r_cmd, demo_dir)
    assert "[1] 235" in r.stdout

    # MAE results
    results_file = "Output/processed_results/mae/mae/MAE_results_v29.tsv"
    r_cmd = """ 
            library(data.table)
            res <- fread("{}")
            print(nrow(res))
            """.format(results_file)
    r = runR(r_cmd, demo_dir)
    assert "335 " in r.stdout


def test_export(demo_dir):
    r = run(["snakemake", "exportCounts", "-j", cores], demo_dir)
    assert "Finished job 0." in r.stderr


def test_all(demo_dir):
    r = run(["snakemake", "-j", cores], demo_dir)
    assert "Finished job 0." in r.stderr
