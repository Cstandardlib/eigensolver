-- xmake.lua

-- Define the project
target("test")

    -- Set the project directory
    -- set_projectdir("./simple")

    -- set_toolchain("gcc")

    set_kind("binary")

    -- Set the language
    -- set_languages("cxx11")
    set_languages("c++17")
    -- set_languages("c++20")

    before_build(function (target)
        -- os.cp("matrix/*", "$(buildir)")
        os.cp("matrix/sparseA.mtx", "$(buildir)")
        os.cp("matrix/a1.mtx", "$(buildir)")
     end)

    -- Add the source files
    -- add_files("test_lobpcg.cpp")

    add_files("tests/simple-test.cpp")

    -- add_files("tests/Eigentest.cpp")
    -- add_files("tests/ortho-test.cpp")
    -- add_files("tests/QRtest.cpp")
    -- add_files("tests/block-test.cpp")
    -- add_files("tests/time-test.cpp")
    -- add_files("tests/matvec-test.cpp")
    -- add_files("tests/sparse-test.cpp")

    -- add_files("lobpcg.cpp")
    -- add_files("matvec.cpp")
    -- add_files("ortho.cpp")

    -- Add the include directories
    add_includedirs("eigen")
    add_includedirs("fast_matrix_market/include")
    add_includedirs(".")

    -- -- Set the build mode and output directory
    -- if is_mode("debug") then
    --     set_symbols("debug")
    --     set_optimize("none")
    --     set_targetdir("$(projectdir)/build/debug")
    -- else
    --     set_optimize("fastest")
    --     set_targetdir("$(projectdir)/build/release")
    -- end