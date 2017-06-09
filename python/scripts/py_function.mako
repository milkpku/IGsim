% for fobj in functions:
  <%
    params = ",\n  ".join(["%s %s" % (t, var) for t, var in fobj])
    inputs = ", ".join([var for t, var in fobj])
    pyargs= ", ".join(["py::arg(\"%s\")" % var for t, var in fobj])
  %>
m.def("${function_name}", []
(
  ${params}
)
{
  return sim::${function_name}(${inputs});
}, __doc_sim_${function_name},
${pyargs}
);

% endfor
