import os
from os.path import isfile, join

def collect_headers(path):
    files = [f for f in os.listdir(path) if isfile(join(path, f))]

    # filename should be "xxxxx.h"
    headers = [f for f in files if f.endswith(".h")]

    return headers

def get_doc(filename, path):
    lines = open(join(path, filename)).readlines()
    non_empty_lines = list(filter(lambda x: x, [l.strip(" \t\n\r") for l in lines]))

    while(non_empty_lines and not non_empty_lines[0].startswith("namespace")):
        non_empty_lines.pop(0)

    docl = filter(
            lambda l: l and l[0] == '/' or l[0] == '*',
            non_empty_lines)

    doc = "\n".join(docl)

    return doc

if __name__=="__main__":

    import sys
    if (len(sys.argv) != 2):
        print("usage: %s <src dir>" % sys.argv[0])

    path = sys.argv[1]
    print("searching docs in %s" % path)
    # get headers
    headers = collect_headers(path)

    # for each header, generate docs
    doc_h = open("py_doc.h", 'w')
    doc_cpp = open("py_doc.cpp", 'w')

    for f in headers:
        print("generating docs for %s" %f)
        doc_h.write("extern const char* __doc_sim_%s;\n" % f[:-2])
        doc_temp = get_doc(f, path)
        doc_cpp.write("const char* __doc_sim_%s = R\"sim_Qu8mg5v7(%s)sim_Qu8mg5v7\";\n" % (f[:-2], doc_temp))

    doc_h.close()
    doc_cpp.close()
