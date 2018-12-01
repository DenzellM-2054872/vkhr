#ifndef VKHR_VULKAN_HAIR_STYLE_HH
#define VKHR_VULKAN_HAIR_STYLE_HH

#include <vkhr/hair_style.hh>
#include <vkhr/vulkan/pipeline.hh>
#include <vkhr/vulkan/drawable.hh>
#include <vkhr/vulkan/depth_map.hh>
#include <vkhr/camera.hh>

#include <vkpp/buffer.hh>
#include <vkpp/command_buffer.hh>
#include <vkpp/descriptor_set.hh>
#include <vkpp/pipeline.hh>

namespace vk = vkpp;

namespace vkhr {
    class Rasterizer;
    namespace vulkan {
        class HairStyle final : public Drawable {
        public:
            HairStyle(const vkhr::HairStyle& hair_style,
                      vkhr::Rasterizer& vulkan_renderer);

            HairStyle() = default;

            void load(const vkhr::HairStyle& hair_style,
                      vkhr::Rasterizer& scene_renderer);

            void draw(vk::CommandBuffer& command_buffer,
                      std::size_t framebuffer) override;

            static void build_pipeline(Pipeline& pipeline_reference,
                                       Rasterizer& vulkan_renderer);

        private:
            vk::IndexBuffer  vertices;
            vk::VertexBuffer positions;
            vk::VertexBuffer tangents;
        };
    }
}

#endif