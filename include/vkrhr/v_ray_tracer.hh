#ifndef VKRHR_VRAYTRACER_HH
#define VKRHR_VRAYTRACER_HH

#include <vkhr/rasterizer/model.hh>
#include <vkhr/rasterizer/hair_style.hh>
#include <vkhr/rasterizer/billboard.hh>
#include <vkhr/rasterizer/linked_list.hh>
#include <vkhr/rasterizer/volume.hh>
#include <vkhr/rasterizer/interface.hh>

#include <vkhr/rasterizer/depth_map.hh>
#include <vkhr/renderer.hh>
#include <vkhr/rasterizer/pipeline.hh>

#include <queue>
#include <vector>
#include <unordered_map>
#include <string>

namespace vk = vkpp;

namespace vkrhr {
    class V_Raytracer final : public vkhr::Renderer {
    public:
        V_Raytracer(vkhr::Window& window, const vkhr::SceneGraph& scene_graph);
        void load(const vkhr::SceneGraph&);
        void draw(const vkhr::SceneGraph&);
    private:
        vk::Instance m_instance;
        vk::Surface m_window_surface;
        vk::PhysicalDevice m_physical_device;
        vk::Device m_device;
        vk::SwapChain m_swap_chain;
        vk::Sampler m_depth_sampler;

        vk::CommandPool m_command_pool;
        vk::DescriptorPool m_descriptor_pool;

        vk::RenderPass m_depth_pass;
        vk::RenderPass m_color_pass;
        vk::RenderPass m_imgui_pass;

        std::vector<vk::Framebuffer> m_framebuffers;
        std::vector<vk::Semaphore> m_image_available, m_render_complete;
        std::vector<vk::Fence> m_command_buffer_finished;

        std::vector<vk::UniformBuffer> m_camera;
        std::vector<vk::UniformBuffer> m_lights;
        std::vector<vk::UniformBuffer> m_params;

        std::vector<vkhr::vulkan::DepthMap> m_shadow_maps;
        std::unordered_map<const vkhr::HairStyle*, vkhr::vulkan::HairStyle> m_hair_styles;
        std::unordered_map<const vkhr::Model*, vkhr::vulkan::Model> m_models;

        vkhr::vulkan::Billboard m_fullscreen_billboard;
        //vkhr::vulkan::LinkedList m_ppll;

        std::vector<vk::QueryPool> m_query_pools;
        std::vector<vk::CommandBuffer> m_command_buffers;

        vkhr::Interface m_imgui;

        void build_render_passes();

        friend class vkhr::vulkan::HairStyle;
        friend class vkhr::vulkan::Model;
        friend class vkhr::vulkan::Volume;
        friend class vkhr::vulkan::Billboard;
        friend class vkhr::vulkan::LinkedList;
        friend class vkhr::vulkan::DepthMap;

        friend class ::vkhr::Interface;
    };
}

#endif
