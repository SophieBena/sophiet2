import numpy as np


def are_close(a, b, abs_tol=1e-8, rel_tol=1e-5):
    return np.allclose(a, b, atol=abs_tol, rtol=rel_tol)


def are_roughly_close(a, b, rel_tol_bounds=[1e-15, 1e-5], abs_tol_bounds=[1e-15, 1e-8]):
    return (not are_close(a, b, rel_tol=rel_tol_bounds[0], abs_tol=abs_tol_bounds[0]) and
            are_close(a, b, rel_tol=rel_tol_bounds[1], abs_tol=abs_tol_bounds[1]))

def all_arrays_are_close(a, b, abs_tol=1e-8, rel_tol=1e-5):
    elems_equal = []
    for arr_a, arr_b in zip(a, b):
        elems_equal.append(np.all(are_close(arr_a, arr_b,
                                            abs_tol=abs_tol, rel_tol=rel_tol)))

    return np.all(elems_equal)

def dict_values_are_close(a, b):
    return a == b
    # assert are_close(np.array(surface_props.values()), np.array(EXPECTED_SURFACE_PROPS.values()),
    #                  rel_tol=1e-6, abs_tol=1e-12)


