import subprocess

num_cells = 2
seed = "NONSENSE"
cil = 60
coa = 24
out_dir = "./output"
snap_period = 1
final_t = 10

subprocess.run(["python",
                "./py_model/main.py",
                "--name", "py_comp_{}_seed={}".format(num_cells, seed),
                "--final_t", final_t, "--snap_period", snap_period,
                "--cil", cil, "--coa", coa,
                "--out_dir", out_dir])
