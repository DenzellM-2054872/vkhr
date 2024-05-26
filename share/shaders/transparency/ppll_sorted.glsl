#ifndef VKHR_PPLL_SORTED_GLSL
#define VKHR_PPLL_SORTED_GLSL


layout(binding = 10, r32ui) uniform uimage2D ppll_heads_sorted;
layout(binding = 11, std430) buffer LinkedList_sorted {
    Node ppll_nodes_sorted[];
};

layout(binding = 12) uniform Config_sorted{ uint ppll_size_sorted; };
layout(binding = 13, std430) buffer LinkedListCounter_sorted {
    uint ppll_counter_sorted;
};

uint ppll_next_node_sorted() {
    uint next_node = atomicAdd(ppll_counter_sorted, 1u);
    if (next_node > ppll_size_sorted)
        return PPLL_NULL_NODE;
    ppll_nodes_sorted[next_node].prev = PPLL_NULL_NODE;
    return next_node;
}

void ppll_node_data_sorted(uint node, vec4 color, float depth) {
    ppll_nodes_sorted[node].color = packUnorm4x8(color);
    ppll_nodes_sorted[node].depth = depth; // don't pack
}

void ppll_node_data_sorted(uint index, Node node) {
    ppll_nodes_sorted[index] = node;
}

Node ppll_node_sorted(uint node) {
    return ppll_nodes_sorted[node];
}

uint ppll_head_node_sorted(ivec2 pixel) {
    return imageLoad(ppll_heads_sorted, pixel).r;
}

void ppll_link_node_sorted(ivec2 pixel, uint next_node) {
    uint prev_node = imageAtomicExchange(ppll_heads_sorted, pixel,
        next_node);
    ppll_nodes_sorted[next_node].prev = prev_node;
}

void ppll_set_head_sorted(ivec2 pixel, uint next_node) {
    imageAtomicExchange(ppll_heads_sorted, pixel, next_node);
}

#endif
