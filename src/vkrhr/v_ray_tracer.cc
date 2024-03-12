#include "vkrhr/v_ray_tracer.hh"

#include <vkpp/append.hh>
#include <vkhr/rasterizer.hh>
namespace vkrhr {
    V_Raytracer::V_Raytracer(vkhr::Window& window, const vkhr::SceneGraph& scene_graph)
    {
        vk::Version target_vulkan_loader{ 1,3 };
        vk::Application application_information{
            "VKHR", { 1, 0, 0 },
            "None", { 0, 0, 0 },
            target_vulkan_loader
        };

        std::vector<vk::Layer> required_layers{
        #ifdef DEBUG
            "VK_LAYER_KHRONOS_validation"
        #endif
        };

        std::vector<vk::Extension> required_extensions{

        #ifdef DEBUG
            "VK_EXT_debug_utils"
        #endif
        };

        vk::append(window.get_vulkan_surface_extensions(),
            required_extensions); // VK_surface_KHR

        m_instance = vk::Instance{
            application_information,
            required_layers,
            required_extensions
        };

        m_window_surface = window.create_vulkan_surface_with(m_instance);

        // Find physical devices that seem most promising of the lot.
        auto score = [&](const vk::PhysicalDevice& physical_device) {
            short gpu_suitable = 2 * physical_device.is_discrete_gpu() +
                physical_device.is_integrated_gpu();
            return physical_device.has_every_queue() * gpu_suitable *
                physical_device.has_present_queue(m_window_surface);
            };

        m_physical_device = m_instance.find_physical_devices_with(score);
        window.append_string(m_physical_device.get_name()); // our GPU.
        m_physical_device.assign_present_queue_indices(m_window_surface);

        
        std::vector<vk::Extension> device_extensions{
            "VK_KHR_acceleration_structure",
            "VK_EXT_descriptor_indexing",
            "VK_KHR_buffer_device_address",
            "VK_KHR_deferred_host_operations",
            "VK_KHR_ray_tracing_pipeline",
            "VK_KHR_spirv_1_4",
            "VK_KHR_shader_float_controls",
            "VK_KHR_ray_query",
            "VK_KHR_swapchain",
        };

        // Just enable every device feature we have right now.
        auto device_features = m_physical_device.get_features();

        m_device = vk::Device{
            m_physical_device,
            required_layers,
            device_extensions,
            device_features
        };

        m_command_pool = vk::CommandPool{ m_device, m_device.get_graphics_queue() };

        auto presentation_mode = vk::SwapChain::mode(window.vsync_requested());

        m_swap_chain = vk::SwapChain{
            m_device,
            m_window_surface,
            m_command_pool,
            {
                VK_FORMAT_B8G8R8A8_UNORM,
                VK_COLOR_SPACE_SRGB_NONLINEAR_KHR
            },
            presentation_mode,
            window.get_extent()
        };

        m_depth_sampler = vk::Sampler{
            m_device, // for sampling depth buffer.
            VK_FILTER_LINEAR,    VK_FILTER_LINEAR,
            VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE,
            VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE
        };

        m_descriptor_pool = vk::DescriptorPool{
            m_device,
            {
                { VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,         64 },
                { VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 64 },
                { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,         64 },
                { VK_DESCRIPTOR_TYPE_INPUT_ATTACHMENT,       64 }
            }
        };

        build_render_passes();

        m_framebuffers = m_swap_chain.create_framebuffers(m_color_pass);

        m_image_available = vk::Semaphore::create(m_device, m_swap_chain.size(), "Image Available Semaphore");
        m_render_complete = vk::Semaphore::create(m_device, m_swap_chain.size(), "Render Complete Semaphore");
        m_command_buffer_finished = vk::Fence::create(m_device, m_swap_chain.size(), "Buffer Finished Fence");

        m_camera = vk::UniformBuffer::create(m_device, sizeof(vkhr::ViewProjection), m_swap_chain.size(), "Camera Matrix Data");
        m_params = vk::UniformBuffer::create(m_device, sizeof(vkhr::Interface::Parameters), m_swap_chain.size(), "Rendering Settings");

        // im hoping i can get away without having to reimplement these...
        // only time shall tell
        //m_ppll = vkhr::vulkan::LinkedList{
        //    *this,
        //    m_swap_chain.get_width(), m_swap_chain.get_height(),
        //    vkhr::vulkan::LinkedList::NodeSize,
        //    vkhr::vulkan::LinkedList::AverageFragmentsPerPixel * m_swap_chain.get_width() *
        //                                                         m_swap_chain.get_height()
        //};

        //fullscreen_billboard = vulkan::Billboard{
        //    swap_chain.get_width(),
        //    swap_chain.get_height(),
        //    *this
        //};

        m_imgui = vkhr::Interface{ m_window_surface.get_glfw_window(), this };

        load(scene_graph);

        m_query_pools = vk::QueryPool::create(m_framebuffers.size(), m_device, VK_QUERY_TYPE_TIMESTAMP, 128);
        m_command_buffers = m_command_pool.allocate(m_framebuffers.size());
    }


    void V_Raytracer::load(const vkhr::SceneGraph& scene_graph)
    {
        m_device.wait_idle();

        m_hair_styles.clear();
        m_models.clear();
        m_shadow_maps.clear();

        for (const auto& model : scene_graph.get_models())
            m_models[&model.second] = vkhr::vulkan::Model{
                model.second, *this
        };

        for (const auto& hair_style : scene_graph.get_hair_styles())
            m_hair_styles[&hair_style.second] = vkhr::vulkan::HairStyle{
                hair_style.second, *this
        };

        m_lights = vk::UniformBuffer::create(m_device, scene_graph.get_light_sources().size() * sizeof(vkhr::LightSource::Buffer),
            m_swap_chain.size(), "Light Source Buffer Data"); // e.g.: position, intensity.

        for (auto& light_source : scene_graph.get_light_sources())
            m_shadow_maps.emplace_back(1024, *this, light_source);

        //build_pipelines();
    }

    void V_Raytracer::draw(const vkhr::SceneGraph&)
    {
    }
    void V_Raytracer::build_render_passes()
    {
        vk::RenderPass::create_modified_color_pass(m_color_pass, m_device, m_swap_chain);
        vk::RenderPass::create_standard_depth_pass(m_depth_pass, m_device);
        vk::RenderPass::create_standard_imgui_pass(m_imgui_pass, m_device, m_swap_chain);
    }
}
