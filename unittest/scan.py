import os
from os.path import isfile, join, dirname, realpath

def collect_headers(path):
    files = [f for f in os.listdir(path) if isfile(join(path, f))]

    # filename should be "<function_name>.h"
    headers = [f[:-2] for f in files if f.endswith(".h")]

    return headers

def collect_unittest(path):
    files = [f for f in os.listdir(path) if isfile(join(path, f))]

    # filename should be "<function_name>.h"
    unittests = [f[:-len("_unittest.cpp")] for f in files if f.endswith(".cpp")]

    return unittests

if __name__=="__main__":

    script_dir = dirname(realpath(__file__))
    header_path = "../include/IGsim"
    unit_path = ''

    headers = collect_headers(join(script_dir, header_path))

    unittests = collect_unittest(join(script_dir, unit_path))

    print("unittest detected:")
    for u in unittests:
        print("-", u)

    print("function need test:")
    for f in (set(headers) - set(unittests)):
        print("-", f)
