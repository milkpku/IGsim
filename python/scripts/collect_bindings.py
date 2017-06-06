import os
from os.path import isfile, join

from mako.template import Template

def collect_functions(path):
    files = [f for f in os.listdir(path) if isfile(join(path, f))]

    # filename should be "py_<function_name>.cpp"
    functions = [f[3:-4] for f in files]

    return functions


if __name__=="__main__":

    path = "../py_sim"

    functions = collect_functions(path)
    sim_temp = Template(filename="py_sim.mako")

    wh = open("py_sim.cpp", 'w')
    wh.write(sim_temp.render(functions = functions))
    wh.close()

    shared_temp = Template(filename="python_shared.mako")
    wh = open("python_shared.cpp", 'w')
    wh.write(shared_temp.render(functions = functions))
    wh.close()
