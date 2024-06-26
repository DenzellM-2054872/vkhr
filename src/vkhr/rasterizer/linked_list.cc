#include <vkhr/rasterizer/linked_list.hh>

#include <vkhr/rasterizer.hh>

namespace vkhr {
    namespace vulkan {
        LinkedList::LinkedList(vkhr::Rasterizer& rasterizer, std::uint32_t width, std::uint32_t height, std::size_t node_size, std::size_t node_count) {
            create(rasterizer, width, height, node_size, node_count);
        }

        void LinkedList::create(vkhr::Rasterizer& rasterizer, std::uint32_t width, std::uint32_t height, std::size_t node_size, std::size_t node_count) {
            heads = vk::DeviceImage {
                rasterizer.device,
                width, height,
                width * height * sizeof(std::uint32_t), VK_FORMAT_R32_UINT,
                VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT
            };

            heads_sorted = vk::DeviceImage{
                rasterizer.device,
                width, height,
                width * height * sizeof(std::uint32_t), VK_FORMAT_R32_UINT,
                VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT
            };

            this->width  = width;
            parameters_buffer.node_count = node_count;
            this->height = height;

            auto command_buffer = rasterizer.command_pool.allocate_and_begin();
            heads.transition(command_buffer,
                             0,
                             VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                             VK_IMAGE_LAYOUT_UNDEFINED,
                             VK_IMAGE_LAYOUT_GENERAL,
                             VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT,
                             VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT);

            heads_sorted.transition(command_buffer,
                0,
                VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                VK_IMAGE_LAYOUT_UNDEFINED,
                VK_IMAGE_LAYOUT_GENERAL,
                VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT,
                VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT);

            command_buffer.end();
            command_buffer.get_queue().submit(command_buffer)
                                      .wait_idle();

            vk::DebugMarker::object_name(rasterizer.device, heads, VK_OBJECT_TYPE_IMAGE, "PPLL Heads", id);
            vk::DebugMarker::object_name(rasterizer.device, heads.get_device_memory(), VK_OBJECT_TYPE_DEVICE_MEMORY,
                                         "PPLL Heads Device Memory", id);

            heads_view = vk::ImageView {
                rasterizer.device,
                heads,
                VK_IMAGE_LAYOUT_GENERAL
            };

            vk::DebugMarker::object_name(rasterizer.device, heads_view, VK_OBJECT_TYPE_IMAGE, "PPLL Heads View", id);

            nodes = vk::StorageBuffer {
                rasterizer.device,
                node_count * node_size
            };

            this->node_size = node_size;

            vk::DebugMarker::object_name(rasterizer.device, nodes, VK_OBJECT_TYPE_BUFFER, "PPLL Nodes", id);
            vk::DebugMarker::object_name(rasterizer.device, nodes.get_device_memory(), VK_OBJECT_TYPE_DEVICE_MEMORY,
                                         "PPLL Nodes Device Memory", id);

            parameters = vk::UniformBuffer {
                rasterizer.device,
                parameters_buffer
            };

            vk::DebugMarker::object_name(rasterizer.device, parameters, VK_OBJECT_TYPE_BUFFER, "PPLL Parameters", id);

            node_counter = vk::StorageBuffer {
                rasterizer.device,
                sizeof(std::uint32_t)
            };

            vk::DebugMarker::object_name(rasterizer.device, node_counter, VK_OBJECT_TYPE_BUFFER, "PPLL Counter", id);

            null_value.uint32[0] = Null;

            vk::DebugMarker::object_name(rasterizer.device, heads_sorted, VK_OBJECT_TYPE_IMAGE, "PPLL sorted Heads", id);
            vk::DebugMarker::object_name(rasterizer.device, heads_sorted.get_device_memory(), VK_OBJECT_TYPE_DEVICE_MEMORY,
                "PPLL Heads Device Memory", id);

            heads_view_sorted = vk::ImageView{
                rasterizer.device,
                heads_sorted,
                VK_IMAGE_LAYOUT_GENERAL
            };

            vk::DebugMarker::object_name(rasterizer.device, heads_view_sorted, VK_OBJECT_TYPE_IMAGE, "PPLL sorted Heads View", id);

            nodes_sorted = vk::StorageBuffer{
                rasterizer.device,
                node_count * node_size
            };

            this->node_size = node_size;

            vk::DebugMarker::object_name(rasterizer.device, nodes_sorted, VK_OBJECT_TYPE_BUFFER, "PPLL sorted Nodes", id);
            vk::DebugMarker::object_name(rasterizer.device, nodes_sorted.get_device_memory(), VK_OBJECT_TYPE_DEVICE_MEMORY,
                "PPLL sorted Nodes Device Memory", id);

            parameters_sorted = vk::UniformBuffer{
                rasterizer.device,
                parameters_buffer
            };

            vk::DebugMarker::object_name(rasterizer.device, parameters_sorted, VK_OBJECT_TYPE_BUFFER, "PPLL sorted Parameters", id);

            node_counter_sorted = vk::StorageBuffer{
                rasterizer.device,
                sizeof(std::uint32_t)
            };

            vk::DebugMarker::object_name(rasterizer.device, node_counter_sorted, VK_OBJECT_TYPE_BUFFER, "PPLL sorted Counter", id);
        }

        void LinkedList::clear(vk::CommandBuffer& command_buffer) {
            heads.transition(command_buffer,
                             VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                             VK_ACCESS_TRANSFER_WRITE_BIT,
                             VK_IMAGE_LAYOUT_GENERAL,
                             VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                             VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                             VK_PIPELINE_STAGE_TRANSFER_BIT);

            command_buffer.clear_color_image(heads, null_value);

            heads.transition(command_buffer,
                             VK_ACCESS_TRANSFER_WRITE_BIT,
                             VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                             VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                             VK_IMAGE_LAYOUT_GENERAL,
                             VK_PIPELINE_STAGE_TRANSFER_BIT,
                             VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT);

            command_buffer.fill_buffer(node_counter,
                                       0,
                                       sizeof(std::uint32_t),
                                       0);

            heads_sorted.transition(command_buffer,
                VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                VK_ACCESS_TRANSFER_WRITE_BIT,
                VK_IMAGE_LAYOUT_GENERAL,
                VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                VK_PIPELINE_STAGE_TRANSFER_BIT);

            command_buffer.clear_color_image(heads_sorted, null_value);

            heads_sorted.transition(command_buffer,
                VK_ACCESS_TRANSFER_WRITE_BIT,
                VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                VK_IMAGE_LAYOUT_GENERAL,
                VK_PIPELINE_STAGE_TRANSFER_BIT,
                VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT);

            command_buffer.fill_buffer(node_counter_sorted,
                0,
                sizeof(std::uint32_t),
                0);
        }

        void LinkedList::resolve(vk::SwapChain& swap_chain, std::uint32_t frame, Pipeline& pipeline, vk::CommandBuffer& command_buffer) {
            swap_chain.get_images()[frame].transition(command_buffer,
                                                      VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT,
                                                      VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                                                      VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL,
                                                      VK_IMAGE_LAYOUT_GENERAL,
                                                      VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT,
                                                      VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT);

            command_buffer.bind_pipeline(pipeline);

            pipeline.descriptor_sets[frame].write(5, heads_view);
            pipeline.descriptor_sets[frame].write(6, nodes);
            pipeline.descriptor_sets[frame].write(7, parameters);
            pipeline.descriptor_sets[frame].write(8, node_counter);
            pipeline.descriptor_sets[frame].write(9, swap_chain.get_general_image_views()[frame]);

            command_buffer.bind_descriptor_set(pipeline.descriptor_sets[frame], pipeline);

            command_buffer.dispatch(std::ceil(width / 8.0), std::ceil(height / 8.0));

            swap_chain.get_images()[frame].transition(command_buffer,
                                                      VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                                                      VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT,
                                                      VK_IMAGE_LAYOUT_GENERAL,
                                                      VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL,
                                                      VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
               
                                       VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT);
        }

        void LinkedList::sort(std::uint32_t frame, Pipeline& pipeline, vk::CommandBuffer& command_buffer) {
            command_buffer.bind_pipeline(pipeline);

            pipeline.descriptor_sets[frame].write(5, heads_view);
            pipeline.descriptor_sets[frame].write(6, nodes);
            pipeline.descriptor_sets[frame].write(7, parameters);
            pipeline.descriptor_sets[frame].write(8, node_counter);

            pipeline.descriptor_sets[frame].write(10, heads_view_sorted);
            pipeline.descriptor_sets[frame].write(11, nodes_sorted);
            pipeline.descriptor_sets[frame].write(12, parameters_sorted);
            pipeline.descriptor_sets[frame].write(13, node_counter_sorted);

            command_buffer.bind_descriptor_set(pipeline.descriptor_sets[frame], pipeline);

            command_buffer.dispatch(std::ceil(width / 8.0), std::ceil(height / 8.0));
        }

        void LinkedList::resolve_Coolchomp(vk::SwapChain& swap_chain, std::uint32_t frame, Pipeline& pipeline, vk::CommandBuffer& command_buffer) {
            swap_chain.get_images()[frame].transition(command_buffer,
                VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT,
                VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL,
                VK_IMAGE_LAYOUT_GENERAL,
                VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT,
                VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT);

            command_buffer.bind_pipeline(pipeline);

            pipeline.descriptor_sets[frame].write(5, heads_view_sorted);
            pipeline.descriptor_sets[frame].write(6, nodes_sorted);
            pipeline.descriptor_sets[frame].write(7, parameters_sorted);
            pipeline.descriptor_sets[frame].write(8, node_counter_sorted);
            pipeline.descriptor_sets[frame].write(9, swap_chain.get_general_image_views()[frame]);

            command_buffer.bind_descriptor_set(pipeline.descriptor_sets[frame], pipeline);

            command_buffer.dispatch(std::ceil(width / 8.0), std::ceil(height / 8.0));

            swap_chain.get_images()[frame].transition(command_buffer,
                VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT,
                VK_IMAGE_LAYOUT_GENERAL,
                VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL,
                VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,

                VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT);
        }
        std::size_t LinkedList::get_width() const {
            return width;
        }

        std::size_t LinkedList::get_node_count() const {
            return parameters_buffer.node_count;
        }

        std::size_t LinkedList::get_node_size() const {
            return node_size;
        }

        std::size_t LinkedList::get_heads_size_in_bytes() const {
            return width * height * sizeof(std::uint32_t);
        }

        vk::UniformBuffer& LinkedList::get_parameters() {
            return parameters;
        }

        std::size_t LinkedList::get_nodes_size_in_bytes() const {
            return parameters_buffer.node_count * node_size;
        }

        std::size_t LinkedList::get_height() const {
            return height;
        }

        void LinkedList::update_resolution(std::size_t width, std::size_t height) {
            this->width  = width;
            this->height = height;
        }

        vk::DeviceImage& LinkedList::get_heads() {
            return heads;
        }

        vk::ImageView& LinkedList::get_heads_view() {
            return heads_view;
        }

        vk::StorageBuffer& LinkedList::get_nodes() {
            return nodes;
        }

        vk::StorageBuffer& LinkedList::get_node_counter() {
            return node_counter;
        }

        void LinkedList::build_pipeline(Pipeline& pipeline, Rasterizer& rasterizer) {
            pipeline = Pipeline { /* In the case we are re-creating the pipeline */ };

            pipeline.shader_stages.emplace_back(rasterizer.device, SHADER("transparency/resolve.comp"));
            vk::DebugMarker::object_name(rasterizer.device, pipeline.shader_stages[0],
                                         VK_OBJECT_TYPE_SHADER_MODULE, "PPLL Resolve");

            pipeline.descriptor_set_layout = vk::DescriptorSet::Layout {
                rasterizer.device,
                {
                    { 5, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE  },
                    { 6, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER },
                    { 7, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER },
                    { 8, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER },
                    { 9, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE  }
                }
            };

            vk::DebugMarker::object_name(rasterizer.device, pipeline.descriptor_set_layout,
                                         VK_OBJECT_TYPE_DESCRIPTOR_SET_LAYOUT, "PPLL Descriptor Set Layout");
            pipeline.descriptor_sets = rasterizer.descriptor_pool.allocate(rasterizer.swap_chain.size(),
                                                                           pipeline.descriptor_set_layout,
                                                                           "PPLL Descriptor Set");

            pipeline.pipeline_layout = vk::Pipeline::Layout {
                rasterizer.device,
                pipeline.descriptor_set_layout
            };

            vk::DebugMarker::object_name(rasterizer.device, pipeline.pipeline_layout,
                                         VK_OBJECT_TYPE_PIPELINE_LAYOUT,
                                         "PPLL Pipeline Layout");

            pipeline.compute_pipeline = vk::ComputePipeline {
                rasterizer.device,
                pipeline.shader_stages[0],
                pipeline.pipeline_layout
            };

            vk::DebugMarker::object_name(rasterizer.device, pipeline.compute_pipeline,
                                         VK_OBJECT_TYPE_PIPELINE, "PPLL sort Pipeline");
        }

        void LinkedList::build_sort_pipeline(Pipeline& pipeline, Rasterizer& rasterizer) {
            pipeline = Pipeline{ /* In the case we are re-creating the pipeline */ };

            pipeline.shader_stages.emplace_back(rasterizer.device, SHADER("transparency/sort_ppll.comp"));
            vk::DebugMarker::object_name(rasterizer.device, pipeline.shader_stages[0],
                VK_OBJECT_TYPE_SHADER_MODULE, "PPLL sort");

            pipeline.descriptor_set_layout = vk::DescriptorSet::Layout{
                rasterizer.device,
                {
                    { 5, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE  },
                    { 6, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER },
                    { 7, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER },
                    { 8, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER },
                    { 10, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE  },
                    { 11, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER },
                    { 12, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER },
                    { 13, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER },
                }
            };

            vk::DebugMarker::object_name(rasterizer.device, pipeline.descriptor_set_layout,
                VK_OBJECT_TYPE_DESCRIPTOR_SET_LAYOUT, "PPLL sort Descriptor Set Layout");
            pipeline.descriptor_sets = rasterizer.descriptor_pool.allocate(rasterizer.swap_chain.size(),
                pipeline.descriptor_set_layout,
                "PPLL sort Descriptor Set");

            pipeline.pipeline_layout = vk::Pipeline::Layout{
                rasterizer.device,
                pipeline.descriptor_set_layout
            };

            vk::DebugMarker::object_name(rasterizer.device, pipeline.pipeline_layout,
                VK_OBJECT_TYPE_PIPELINE_LAYOUT,
                "PPLL sort Pipeline Layout");

            pipeline.compute_pipeline = vk::ComputePipeline{
                rasterizer.device,
                pipeline.shader_stages[0],
                pipeline.pipeline_layout
            };

            vk::DebugMarker::object_name(rasterizer.device, pipeline.compute_pipeline,
                VK_OBJECT_TYPE_PIPELINE, "PPLL sort Pipeline");
        }

    void LinkedList::build_pipeline_Coolchomp(Pipeline& pipeline, Rasterizer& rasterizer) {
        pipeline = Pipeline{ /* In the case we are re-creating the pipeline */ };

        pipeline.shader_stages.emplace_back(rasterizer.device, SHADER("transparency/resolve_Coolchomp.comp"));
        vk::DebugMarker::object_name(rasterizer.device, pipeline.shader_stages[0],
            VK_OBJECT_TYPE_SHADER_MODULE, "PPLL Resolve frost");

        pipeline.descriptor_set_layout = vk::DescriptorSet::Layout{
            rasterizer.device,
            {
                    { 5, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE  },
                    { 6, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER },
                    { 7, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER },
                    { 8, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER },
                    { 9, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE  }
            }

        };
        vk::DebugMarker::object_name(rasterizer.device, pipeline.descriptor_set_layout,
            VK_OBJECT_TYPE_DESCRIPTOR_SET_LAYOUT, "PPLL Descriptor Set Layout frost");
        pipeline.descriptor_sets = rasterizer.descriptor_pool.allocate(rasterizer.swap_chain.size(),
            pipeline.descriptor_set_layout,
            "PPLL Descriptor Set frost");


        pipeline.pipeline_layout = vk::Pipeline::Layout{
            rasterizer.device,
            pipeline.descriptor_set_layout
        };

        vk::DebugMarker::object_name(rasterizer.device, pipeline.pipeline_layout,
            VK_OBJECT_TYPE_PIPELINE_LAYOUT,
            "PPLL Pipeline Layout frost");

        pipeline.compute_pipeline = vk::ComputePipeline{
            rasterizer.device,
            pipeline.shader_stages[0],
            pipeline.pipeline_layout
        };

        vk::DebugMarker::object_name(rasterizer.device, pipeline.compute_pipeline,
            VK_OBJECT_TYPE_PIPELINE, "PPLL Pipeline frost");
    }

        int LinkedList::id { 0 };

    }
 }
