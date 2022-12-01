# a dictionnary containing variable names as keys

# each entry is a dictionnary containing the key "deps", which is a
# list of the variable this variable depends on. example:

# dep_dic["a_variable"]["deps"] gives the list of all the variables
# needed to compute "a_variable"
dep_dic = {"grad_v": {"deps": ["v"]},
           "div_v": {"deps": ["v"]},
           "lap_v": {"deps": ["v"]},
           "vJ": {"deps": ["v", "J"]},
           "div_vJ": {"deps": ["v",
                               "J"]},
           "grad_div_v": {"deps": ["div_v"]},
           "grad_ln_rho_traceless_grad_v": {"deps": ["traceless_grad_v",
                                                     "grad_ln_rho"]},
           "grad_ln_rho": {"deps": ["ln_rho"]},
           "ln_rho": {"deps": ["rho"]},
           "grad_rho": {"deps": ["rho"]},
           "div_press": {"deps": ["press"]},
           "lap_rho": {"deps": ["rho"]},
           "press": {"deps": ["rho", "grad_rho", "lap_rho", "T"]},
           "traceless_grad_v_grad_v": {"deps": ["traceless_grad_v",
                                                "grad_v"]},
           "grad_ln_rho_grad_T": {"deps":["grad_ln_rho",
                                          "grad_T"]},
           "lap_T": {"deps":["T"]},
           "v_grad_T": {"deps":["v",
                                "grad_T"]},
           "grad_T": {"deps":["T"]},
           "v_grad_ln_rho": {"deps":["v", "grad_ln_rho"]}}

# returns a list of all the dependencies required to compute "key" in
# the dependencies dictionnary "depdic". You can keep the doublons if
# you prefer with "doublons" optional variable
def get_dependencies(depdic, key, doublons = False):
    all_dependencies = []
    to_check_dep = [key]
    checked_dep = []
    while len(to_check_dep) != 0:
        for depkey in to_check_dep:
            if depkey in depdic and not(depkey in checked_dep):
                deplist = depdic[depkey]["deps"]
                all_dependencies += deplist
                if not(doublons):
                    all_dependencies = list(set(all_dependencies))
                checked_dep.append(depkey)
                for d in deplist:
                    if not(d in checked_dep):
                        to_check_dep.append(d)
                to_check_dep.remove(depkey)
            else:
                checked_dep.append(depkey)
                to_check_dep.remove(depkey)
    return all_dependencies

# ////////////////////////////////////////////
# TESTS
# ////////////////////////////////////////////

a = get_dependencies(dep_dic, "grad_ln_rho")
a = get_dependencies(dep_dic, "v_grad_ln_rho")
a = get_dependencies(dep_dic, "press")
# print(a)

deps_list = []
for d in dep_dic:
    all_deps = get_dependencies(dep_dic, d)
    deps_list.append([d, len(all_deps)])
    print(f"var: {d}\ndeps: {all_deps}\n---")
sorted_deps_list = sorted(deps_list, key = lambda x: x[1],
                          reverse=True)
print(sorted_deps_list)
    
l = [['div_press', 5],
     ['grad_ln_rho_grad_T', 5],
     ['grad_ln_rho_traceless_grad_v', 4],
     ['press', 4],
     ['v_grad_ln_rho', 4],
     ['traceless_grad_v_grad_v', 3],
     ['v_grad_T', 3],
     ['div_vJ', 2],
     ['grad_div_v', 2],
     ['grad_ln_rho', 2],
     ['grad_v', 1],
     ['div_v', 1],
     ['lap_v', 1],
     ['ln_rho', 1],
     ['grad_rho', 1],
     ['lap_rho', 1],
     ['lap_T', 1],
     ['grad_T', 1]]
