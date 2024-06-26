#ifndef VKHR_VULKAN_HAIR_STYLE_HH
#define VKHR_VULKAN_HAIR_STYLE_HH

#include <vkhr/scene_graph/hair_style.hh>

#include <vkhr/rasterizer/pipeline.hh>
#include <vkhr/rasterizer/volume.hh>
#include <vkhr/rasterizer/drawable.hh>

#include <vkhr/scene_graph/camera.hh>

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
            //void gen_DOM    (Pipeline& DOM_pipeline,            vk::DescriptorSet& descriptor_set, vk::CommandBuffer& command_buffer);
            void voxelize   (Pipeline& voxelization_pipeline,   vk::DescriptorSet& descriptor_set, vk::CommandBuffer& command_buffer);
            void draw_volume(Pipeline& volume_pipeline,         vk::DescriptorSet& descriptor_set, vk::CommandBuffer& command_buffer);
            //void voxelize_Coolchomp(Pipeline& voxelization_pipeline, vk::DescriptorSet& descriptor_set, vk::CommandBuffer& command_buffer);

            void draw(Pipeline& vulkan_strand_rasterizer_pipeline,
                      vk::DescriptorSet& descriptor_set,
                      vk::CommandBuffer& command_buffer) override;
            static void build_Coolchomp_pipeline(Pipeline& pipeline, Rasterizer& vulkan_renderer);
            //static void dom_pipeline(Pipeline& pipeline_reference, Rasterizer& vulkan_renderer, int dom_layers);

            static void build_pipeline(Pipeline& pipeline_reference, Rasterizer& vulkan_renderer);
            static void depth_pipeline(Pipeline& pipeline_reference, Rasterizer& vulkan_renderer);
            static void voxel_pipeline(Pipeline& pipeline_reference, Rasterizer& vulkan_renderer);
            static void voxel_pipeline_Coolchomp(Pipeline& pipeline_reference, Rasterizer& vulkan_renderer);
            
            void update_parameters();

            std::size_t get_geometry_size() const;
            std::size_t get_volume_size()   const;

            struct Parameters {
                AABB volume_bounds;
                glm::vec3 volume_resolution;
                float strand_radius;
                glm::vec3 hair_color;
                float hair_opacity;
                float hair_shininess;
                float strand_ratio;
                int vertex_count;
            } parameters;

            void reduce(float ratio);

        private:

            uint32_t dom_width;
            uint32_t dom_height;

            vk::IndexBuffer  segments;
            vk::VertexBuffer vertices;
            vk::VertexBuffer tangents;
            vk::VertexBuffer thickness;

/*          vk::ImageView transmitance_view;
            vk::DeviceImage transmitance_volume;
            vk::Sampler transmitance_sampler;

            vk::ImageView variance_view;
            vk::DeviceImage variance_volume;
            vk::Sampler variance_sampler;
*/

            vk::ImageView density_view;
            vk::DeviceImage density_volume;
            vk::Sampler density_sampler;

            //std::vector<vk::ImageView> dom_views;
            //std::vector<vk::DeviceImage> dom_volumes;
            //std::vector<vk::Sampler> dom_samplers;

            vk::ImageView tangent_view;
            vk::DeviceImage tangent_volume;
            vk::Sampler tangent_sampler;

            vk::UniformBuffer parameter_buffer;

            Volume volume;

            std::size_t segments_per_strand;

            friend class Volume;

            static int id;
        };
    }
}

#endif