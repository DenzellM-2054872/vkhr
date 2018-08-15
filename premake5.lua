name = "vkhr"
workspace (name)
    language "C++"
    location "build"
    warnings "Extra"
    cppdialect "C++17"
    configurations { "Debug",
                     "Release" }

    filter { "action:vs*" }
        startproject "vkhr"
        architecture "x86_64"
        systemversion "10.0.17134.0"

    filter "configurations:Debug"
        floatingpoint "Strict"
        defines "DEBUG"
        optimize "Off"
        symbols "On"

    filter "configurations:Release"
        floatingpoint "Fast"
        defines "RELEASE"
        optimize "Speed"
        symbols "Off"

SDK = "$(VULKAN_SDK)"

project (name)
    targetdir "bin"
    kind "ConsoleApp"

    includedirs "include"
    files { "include/**.hh" }
    files { "src/"..name.."/**.cc" }
    files "src/main.cc"

    files { "share/scenes/**.vkhr" }
    files { "share/shader/**.glsl",
            "share/shader/**.vert",
            "share/shader/**.geom",
            "share/shader/**.tesc",
            "share/shader/**.tese",
            "share/shader/**.frag",
            "share/shader/**.comp" }

    vpaths {
        ["src/*"] = "src/**.cc",
        ["include/*"] = "include/**.hh",
        ["scenes/*"]  = { "share/scenes/**.vkhr" }
        ["shaders/*"] = { "share/shader/**.glsl",
                          "share/shader/**.vert",
                          "share/shader/**.geom",
                          "share/shader/**.tesc",
                          "share/shader/**.tese",
                          "share/shader/**.frag",
                          "share/shader/**.comp" }
    }

    -- For header-only libraries.
    includedirs "foreign/include"
    includedirs "foreign/glm"

    filter { "system:windows", "action:gmake" }
        links { "glfw3", "vulkan" }
    filter { "system:windows", "action:vs*" }
        links { SDK.."/lib/vulkan-1.lib" }
        includedirs { SDK.."/include" }
    filter "system:linux or bsd or solaris"
        links { "glfw3", "vulkan" }