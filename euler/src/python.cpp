#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <euler/euler.h>

namespace py = pybind11;


PYBIND11_MAKE_OPAQUE(std::vector<var>);
PYBIND11_MAKE_OPAQUE(std::vector<inputVar>);

PYBIND11_MODULE(euler, m) {
    m.doc() = "Euler equations solver";

    py::class_<constants>(m, "constants")
        .def(py::init<const real>(), py::arg("gamma"))
        .def(py::init<>())
        .def_readwrite("gamma", &constants::gamma)
        ;

    py::class_<vec2>(m, "vec2")
        .def(py::init<const real, const real>())
        .def(py::init<>())
        .def_readwrite("x", &vec2::x)
        .def_readwrite("y", &vec2::y)
        ;

    py::class_<inputVar>(m, "var")
        .def(py::init<const real, const real, const real, const real>(),
            py::arg("rho"), py::arg("u"), py::arg("v"), py::arg("p")
        )
        .def(py::init<>())
        .def_readwrite("rho", &inputVar::rho)
        .def_readwrite("u", &inputVar::u)
        .def_readwrite("v", &inputVar::v)
        .def_readwrite("p", &inputVar::p)
        .def("__getitem__", [](inputVar& v, uint index) { return &v[index];})
        .def("__len__", &inputVar::size)
        ;

    py::class_<var>(m, "rhovar")
        .def(py::init<const real, const real, const real, const real>())
        .def(py::init<>())
        .def_readwrite("rho", &var::rho)
        .def_readwrite("rhou", &var::rhou)
        .def_readwrite("rhov", &var::rhov)
        .def_readwrite("rhoe", &var::rhoe)
        .def("__getitem__", [](var& v, uint index) { return &v[index];})
        .def("__len__", &var::size)
        ;

    py::bind_vector<std::vector<var>>(m, "vectorRhoVar");
    py::bind_vector<std::vector<inputVar>>(m, "vectorVar");
    
    py::class_<mesh>(m, "mesh")
        .def(py::init<>())
        .def(py::init<const std::string, const bool>(), py::arg("filename"), py::arg("verb")=true)
        .def("size", &mesh::size)
        .def("__len__", &mesh::size)
        .def("get_xy", &mesh::get_xy)
        ;

    py::class_<solver>(m, "solver")
        .def("set_mesh", &solver::set_mesh)
        .def("set_cfl", &solver::set_cfl)
        .def("set_max_time", &solver::set_max_time)
        .def("set_print_interval", &solver::set_print_interval)
        .def("set_max_steps", &solver::set_max_steps)
        .def("set_mesh", &solver::set_mesh)
        .def("set_order", &solver::set_order)
        .def("set_coefficient_save_field", &solver::set_coefficient_save_field)
        .def("set_tolerance", &solver::set_tolerance)
        .def("set_smoother", &solver::set_smoother,
            py::arg("kind"), py::arg("coeff") = -1., py::arg("iterations") = 2
        )
        .def("set_limiter", &solver::set_limiter)
        .def("set_bc", &solver::set_bc,
            py::arg("field"), py::arg("kind"), py::arg("value") = inputVar()
        )
        .def("simulate", &solver::simulate, py::arg("verb") = true)
        .def("writeVtk", &solver::writeVtk)
        .def("writeField", &solver::writeField)
        .def("writeResiduals", &solver::writeResiduals)
        .def("writeSavedCoefficients", &solver::writeSavedCoefficients)
        .def("get_residuals", &solver::get_residuals)
        .def("getField", &solver::getField)
        .def("forceOnField", &solver::forceOnField)
        .def("computeCoefficients", &solver::computeCoefficients)
        .def("get_solution", &solver::get_solution)
        .def("get_limiters", &solver::get_limiters)
        .def("log", &solver::log)
        ;
    
    py::class_<transientEulerSolver, solver>(m, "transientEulerSolver")
        .def(py::init<const mesh&, const std::vector<inputVar>&, const constants>(),
            py::arg("mesh"), py::arg("init"), py::arg("constants")
        );

    py::class_<transientRk5Solver, solver>(m, "transientRk5Solver")
        .def(py::init<const mesh&, const std::vector<inputVar>&, const constants>(),
            py::arg("mesh"), py::arg("init"), py::arg("constants")
        );

    py::class_<steadyRk5Solver, solver>(m, "steadyRk5Solver")
        .def(py::init<const mesh&, const std::vector<inputVar>&, const constants>(),
            py::arg("mesh"), py::arg("init"), py::arg("constants")
        );
}





