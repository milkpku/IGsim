from os.path import isfile, basename, dirname, realpath, join
import re

from mako.template import Template

# return list of function_obj for given function file <function>.h
# function_obj = [(type, var_name) ...]
def parse(filename):
    s = open(filename).read()

    # find function definition void <function>(type vars, ...);
    pattern = re.compile(r'void [^(]*\(([^)]*)\)')
    items = pattern.findall(s)

    functions = []
    # for each overloading of function
    for i in items:
        fobj = []
        # split variables by ,
        variables = list(map(lambda x: x.strip(), i.split(",")))
        for v in variables:
            # get type and variable name
            atoms = v.split(" ")
            fobj.append((
                " ".join(atoms[:-1]),
                atoms[-1]))
        functions.append(fobj)

    return functions

if __name__=="__main__":
    import sys

    if (len(sys.argv) != 2):
        print("usage: %s <header_file>" % sys.argv[0])
        exit()

    if (not isfile(sys.argv[1])):
        print("unable to find file: %s" % sys.argv[1])
        exit()

    # parse function file
    function_name = basename(sys.argv[1])[:-2] # eliminate ".h" suffix
    functions = parse(sys.argv[1])

    print("detect %d function overloading for %s" % (len(functions), function_name))
    for fobj in functions:
        print("void %s(\n\t%s)" % (
            function_name,
            ",\n\t".join(["%s %s" % (t, v) for t, v in fobj])))

    script_dir = dirname(realpath(__file__))
    func_temp = Template(filename= join(script_dir ,"py_function.mako"))

    wh = open("py_%s.cpp" % function_name, 'w')
    wh.write(func_temp.render(
        function_name = function_name,
        functions = functions))
    wh.close()
