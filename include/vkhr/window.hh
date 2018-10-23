#ifndef VKHR_WINDOW_HH
#define VKHR_WINDOW_HH

#include <vkhr/image.hh>

#include <vkpp/extension.hh>

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>

#include <string>
#include <vector>

namespace vkhr {
    class Window final {
    public:
        Window(const int width, const int height, const std::string& title,
               const Image& icon, const bool fullscreen = false,
               const bool vsync = false); // Will be Vulkan set.
        ~Window() noexcept;

        bool is_open() const;
        void close(); // kill

        bool is_fullscreen() const;
        bool request_vsync() const;

        void toggle_fullscreen();
        void hide(); void show();

        int get_width()  const;
        float get_aspect_ratio() const;
        int get_height() const;

        int get_refresh_rate() const;

        std::vector<vkpp::Extension> get_surface_extensions() const;
        VkSurfaceKHR create_surface(VkInstance instance) const;

        void resize(const int width, const int height);

        GLFWwindow* get_handle();

        void set_icon(const Image& icon);

        void change_title(const std::string& title);
        void append_string(const std::string& text);

        void poll_events(); // Called once per frame.

        float get_vertical_dpi()   const;
        float get_horizontal_dpi() const;

        void set_time(const float time);
        float get_current_time() const;
        float delta_time() const;

    private:
        int width, height;
        std::string title;
        bool fullscreen { false },
             vsync;

        GLFWwindow* handle { nullptr };

        GLFWmonitor* monitor;
        int window_x, window_y;
        int monitor_width, monitor_height;
        int monitor_refresh_rate;
        float horizontal_dpi,
              vertical_dpi;

        std::string append;

        std::size_t frames { 0 };

        float frame_time,
              last_frame_time,
              fps_update;
    };
}

#endif