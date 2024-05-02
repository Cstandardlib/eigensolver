-- xmake.lua

-- Define the project
target("test")

    set_kind("binary")

    -- Set the language
    set_languages("cxx11")

    -- Add the source files
    add_files("B-ortho-test.cpp")

    -- Add the include directories
    add_includedirs("../eigen")

    -- -- Set the build mode and output directory
    -- if is_mode("debug") then
    --     set_symbols("debug")
    --     set_optimize("none")
    --     set_targetdir("$(projectdir)/build/debug")
    -- else
    --     set_optimize("fastest")
    --     set_targetdir("$(projectdir)/build/release")
    -- end