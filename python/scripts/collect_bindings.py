import os
from os.path import isfile, join, dirname, realpath

from mako.template import Template

def collect_functions(path):
    files = [f for f in os.listdir(path) if isfile(join(path, f))]

    # filename should be "py_<function_name>.cpp"
    functions = [f[3:-4] for f in files]

    return functions

def collect_headers(path):
    files = [f for f in os.listdir(path) if isfile(join(path, f))]

    # filename should be "<function_name>.h"
    headers = [f[:-2] for f in files if f.endswith(".h")]

    return headers

if __name__=="__main__":

    script_dir = dirname(realpath(__file__))
    py_path = "../py_sim"
    header_path = "../../include/IGsim"

    functions = collect_functions(join(script_dir, py_path))
    headers = collect_headers(join(script_dir, header_path))

    print("function detected:")
    for f in functions:
        print("-", f)

    print("function need implement:")
    for f in (set(headers) - set(functions)):
        print("-", f)

    sim_temp = Template(filename=join(script_dir, "py_sim.mako"))
    wh = open("py_sim.cpp", 'w')
    wh.write(sim_temp.render(functions = functions))
    wh.close()

    shared_temp = Template(filename=join(script_dir,"python_shared.mako"))
    wh = open("python_shared.cpp", 'w')
    wh.write(shared_temp.render(functions = functions))
    wh.close()
