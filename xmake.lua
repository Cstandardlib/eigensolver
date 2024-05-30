-- xmake.lua

-- Define the project
target("main")
    set_kind("binary")

    set_languages("c++17")

    set_optimize("fastest")

    -- add sparse matrix files to build directories
    before_build(function (target)
        -- os.cp("matrix/*", "$(buildir)")
        print("making build directory: ./build")
        os.mkdir("build")
        print("Copying sparse matrix files to build directory")
        os.cp("matrix/sparseA.mtx", "$(buildir)")
        os.cp("matrix/a1.mtx", "$(buildir)")
        os.cp("matrix/b1.mtx", "$(buildir)")
        os.cp("matrix/Si5H12.mtx", "$(buildir)")
        os.cp("matrix/Si2.mtx", "$(buildir)")
        os.cp("matrix/Na5.mtx", "$(buildir)")
        os.cp("matrix/Ga3As3H12.mtx", "$(buildir)")
     end)

    -- Add the test files
    add_files("test_lobpcg.cpp")

    -- Add the source files
    add_files("lobpcg.cpp")
    add_files("matvec.cpp")
    add_files("sparseMv.cpp")
    add_files("ortho.cpp")
    add_files("utils.cpp")

    -- Add the include directories
    add_includedirs("eigen")
    add_includedirs("fast_matrix_market/include")
    add_includedirs(".")

target("test")

    -- Set the project directory
    -- set_projectdir("./simple")

    -- set_toolchain("gcc")

    set_kind("binary")

    -- Set the language
    -- set_languages("cxx11")
    set_languages("c++17")

    -- add sparse matrix files to build directories
    before_build(function (target)
        -- os.cp("matrix/*", "$(buildir)")
        os.cp("matrix/sparseA.mtx", "$(buildir)")
        os.cp("matrix/a1.mtx", "$(buildir)")
        os.cp("matrix/b1.mtx", "$(buildir)")
     end)

    -- Add the test files
    -- add_files("test_lobpcg.cpp")

    -- add_files("tests/simple-test.cpp")

    -- add_files("tests/eigen-test.cpp")
    add_files("tests/ortho-test.cpp")
    -- add_files("tests/QRtest.cpp")
    -- add_files("tests/block-test.cpp")
    -- add_files("tests/time-test.cpp")
    -- add_files("tests/matvec-test.cpp")
    -- add_files("tests/sparse-test.cpp")

    -- Add the source files
    -- add_files("lobpcg.cpp")
    -- add_files("matvec.cpp")
    add_files("ortho.cpp")
    -- add_files("utils.cpp")

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