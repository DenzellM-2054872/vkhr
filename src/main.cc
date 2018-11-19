#include <vkhr/vkhr.hh>
#include <vkpp/vkpp.hh>

namespace vk = vkpp;

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#define GLM_ENABLE_EXPERIMENTAL

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

int main(int argc, char** argv) {
    vkhr::ArgParser argp { vkhr::arguments };
    auto scene_file = argp.parse(argc, argv);

    vk::Version target_vulkan_loader { 1,1 };
    vk::Application application_information {
        "VKHR", { 1, 0, 0 },
        "None", { 0, 0, 0 },
        target_vulkan_loader
    };

    std::vector<vk::Layer> required_layers {
        "VK_LAYER_LUNARG_standard_validation"
    };

    std::vector<vk::Extension> required_extensions {
        "VK_EXT_debug_utils"
    };

    const vkhr::Image vulkan_icon { IMAGE("vulkan.icon") };
    vkhr::Window window { 1280, 720, "VKHR", vulkan_icon };

    vkhr::InputMap input_map { window };

    input_map.bind("quit", vkhr::Input::Key::Escape);
    input_map.bind("grab", vkhr::Input::MouseButton::Left);
    input_map.bind("screenshot", vkhr::Input::Key::S);
    input_map.bind("hide", vkhr::Input::Key::H);

    vk::append(window.get_vulkan_surface_extensions(),
               required_extensions); // VK_surface_KHR

    vk::Instance instance {
        application_information,
        required_layers,
        required_extensions
    };

    auto window_surface = window.create_vulkan_surface(instance);

    // Find physical devices that seem most promising of the lot.
    auto score = [&](const vk::PhysicalDevice& physical_device) {
        short gpu_suitable = physical_device.is_discrete_gpu()*2+
                             physical_device.is_integrated_gpu();
        return physical_device.has_every_queue() * gpu_suitable *
               physical_device.has_present_queue(window_surface);
    };

    auto physical_device = instance.find_physical_devices(score);

    window.append_string(physical_device.get_name());

    physical_device.assign_present_queue_indices(window_surface);

    std::vector<vk::Extension> device_extensions {
        "VK_KHR_swapchain"
    };

    // Just enable every device feature we have right now.
    auto device_features = physical_device.get_features();

    vk::Device device {
        physical_device,
        required_layers,
        device_extensions,
        device_features
    };

    vk::CommandPool command_pool { device, device.get_graphics_queue() };

    vk::SwapChain swap_chain {
        device,
        window_surface,
        {
            VK_FORMAT_B8G8R8A8_UNORM,
            VK_COLOR_SPACE_SRGB_NONLINEAR_KHR
        },
        vk::SwapChain::PresentationMode::Fifo,
        window.get_extent()
    };

    std::vector<vk::RenderPass::Attachment> attachments {
        {
            swap_chain.get_color_attachment_format(),
            swap_chain.get_khr_presentation_layout()
        },
        {
            swap_chain.get_depth_attachment_format(),
            swap_chain.get_depth_attachment_layout()
        }
    };

    std::vector<vk::RenderPass::Subpass> subpasses {
        {
            { 0, swap_chain.get_color_attachment_layout() },
            { 1, swap_chain.get_depth_attachment_layout() }
        }
    };

    std::vector<vk::RenderPass::Dependency> dependencies {
        {
            VK_SUBPASS_EXTERNAL,
            0,
            VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT,
            0,
            VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT,
            VK_ACCESS_COLOR_ATTACHMENT_READ_BIT |
            VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT
        }
    };

    vk::RenderPass render_pass {
        device,
        attachments,
        subpasses,
        dependencies
    };

    vkhr::HairStyle hair_style { STYLE("ponytail.hair") };

    vk::VertexBuffer vertex_buffer {
        device,
        command_pool,
        hair_style.get_vertices(),
        0,
        {
            { 0, VK_FORMAT_R32G32B32_SFLOAT }
        }
    };

    hair_style.generate_tangents();

    vk::VertexBuffer tangent_buffer {
        device,
        command_pool,
        hair_style.get_tangents(),
        1,
        {
            { 1, VK_FORMAT_R32G32B32_SFLOAT }
        }
    };

    hair_style.generate_indices();

    vk::IndexBuffer index_buffer {
        device,
        command_pool,
        hair_style.get_indices()
    };

    vk::GraphicsPipeline::FixedFunction fixed_functions;

    fixed_functions.add_vertex_input(vertex_buffer);
    fixed_functions.add_vertex_input(tangent_buffer);

    fixed_functions.set_topology(VK_PRIMITIVE_TOPOLOGY_LINE_LIST);

    fixed_functions.set_scissor({ 0, 0, swap_chain.get_extent() });
    fixed_functions.set_viewport({ 0.0, 0.0,
                                   static_cast<float>(swap_chain.get_width()),
                                   static_cast<float>(swap_chain.get_height()),
                                   0.0, 1.0 });

    fixed_functions.set_line_width(1.0);
    fixed_functions.enable_depth_test();

    fixed_functions.enable_alpha_blending_for(0);

    std::vector<vk::ShaderModule> shading_stages;

    shading_stages.emplace_back(device, SHADER("kajiya-kay.vert"));
    shading_stages.emplace_back(device, SHADER("kajiya-kay.frag"));

    struct Transform {
        glm::mat4 model      { 1.0f };
        glm::mat4 view       { 1.0f };
        glm::mat4 projection { 1.0f };
    } mvp;

    std::vector<vk::UniformBuffer> uniform_buffers;

    for (std::size_t i { 0 } ; i < swap_chain.size(); ++i)
        uniform_buffers.emplace_back(device, mvp,
                                     sizeof(mvp));

    vk::DescriptorSet::Layout descriptor_layout {
        device,
        {
            { 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER }
        }
    };

    vk::DescriptorPool descriptor_pool {
        device,
        {
            { VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, swap_chain.size() },
            { VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1 } // ImGUI
        }
    };

    auto descriptor_sets = descriptor_pool.allocate(swap_chain.size(), descriptor_layout);

    for (std::size_t i { 0 }; i < descriptor_sets.size(); ++i)
        descriptor_sets[i].write(0, uniform_buffers[i]);

    vk::Pipeline::Layout layout {
        device,
        descriptor_layout
    };

    vk::GraphicsPipeline graphics_pipeline {
        device,
        shading_stages,
        fixed_functions,
        layout,
        render_pass
    };

    auto framebuffers = swap_chain.create_framebuffers(render_pass);

    vk::Semaphore image_available { device };
    vk::Fence command_buffer_done { device };
    vk::Semaphore render_complete { device };

    vkhr::Interface imgui {
        window,
        instance,
        device,
        descriptor_pool,
        render_pass,
        command_pool
    };

    auto command_buffers = command_pool.allocate(framebuffers.size());

    vkhr::Camera camera { glm::radians(45.0f), swap_chain.get_width(),
                                               swap_chain.get_height() };
    camera.look_at({ 0.000f, 60.0f, 0.000f }, { 200.0f, 35.0f, 200.0f });

    mvp.projection  = camera.get_projection_matrix();

    glm::vec2 previous_mouse_position { input_map.get_mouse_position() };

    vkhr::Raytracer raytracer { camera, hair_style };

    window.show();

    while (window.is_open()) {
        if (input_map.just_pressed("quit")) {
            window.close();
        } else if (input_map.just_pressed("hide")) {
            imgui.toggle_visibility();
        } else if (input_map.just_pressed("screenshot")) {
            raytracer.draw(camera);
        }

        glm::vec2 cursor_delta { 0.0f, 0.0f };
        if (input_map.just_released("grab")) {
            input_map.unlock_cursor();
        } else if (input_map.just_pressed("grab")) {
            input_map.freeze_cursor();
            previous_mouse_position  = input_map.get_mouse_position();
        } else if (input_map.pressed("grab") && !imgui.want_focus()) {
            glm::vec2 mouse_position = input_map.get_mouse_position();
            cursor_delta = (mouse_position - previous_mouse_position);
            previous_mouse_position = mouse_position; // single frame.
            camera.arcball_by(cursor_delta*window.delta_time()*0.32f);
        }

        imgui.update();

        auto next_image = swap_chain.acquire_next_image(image_available);

        command_buffer_done.wait_and_reset();

        for (std::size_t i = 0; i < framebuffers.size(); ++i) {
            command_buffers[i].begin();
            command_buffers[i].begin_render_pass(render_pass, framebuffers[i],
                                                 { 1.0f, 1.0f, 1.0f, 1.0f });
            command_buffers[i].bind_pipeline(graphics_pipeline);
            command_buffers[i].bind_descriptor_set(descriptor_sets[i],
                                                   graphics_pipeline);
            command_buffers[i].bind_index_buffer(index_buffer);
            command_buffers[i].bind_vertex_buffer(tangent_buffer);
            command_buffers[i].bind_vertex_buffer(vertex_buffer);
            command_buffers[i].draw_indexed(index_buffer.count());
            imgui.render(command_buffers[i]); // Re-uploads data.
            command_buffers[i].end_render_pass();
            command_buffers[i].end();
        }

        mvp.view   =   camera.get_view_matrix();
        uniform_buffers[next_image].update(mvp);

        device.get_graphics_queue().submit(command_buffers[next_image], image_available,
                                           VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT,
                                           render_complete, command_buffer_done);

        device.get_present_queue().present(swap_chain, next_image, render_complete);

        window.poll_events();
    }

    return 0;
}
