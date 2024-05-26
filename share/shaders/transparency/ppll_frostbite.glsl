#ifndef VKHR_PPLL_GLSL
#define VKHR_PPLL_GLSL

#define PPLL_NULL_NODE 0xffffffff
#define PPLL_MAXIMUM_DEPTH 3.e+38

struct Node {
    uint position;
    uint tangent;
    uint hair_color;
    float thickness;
    float depth;
    uint prev;
};

layout(binding = 5, r32ui) uniform uimage2D ppll_heads;
layout(binding = 6, std430) buffer LinkedList {
    Node ppll_nodes[];
};

layout(binding = 7) uniform Config{ uint ppll_size; };
layout(binding = 8, std430) buffer LinkedListCounter {
    uint ppll_counter;
};

uint ppll_next_node() {
    uint next_node = atomicAdd(ppll_counter, 1u);
    if (next_node > ppll_size)
        return PPLL_NULL_NODE;
    ppll_nodes[next_node].prev = PPLL_NULL_NODE;
    return next_node;
}

void ppll_node_data(uint node, vec4 position, vec3 tangent, float strand_width, vec3 hair_color, float hair_alpha, float thickness, float depth) {
    ppll_nodes[node].position = packSnorm4x8(position);
    ppll_nodes[node].tangent = packSnorm4x8(vec4(tangent, strand_width));
    ppll_nodes[node].hair_color = packSnorm4x8(vec4(hair_color, hair_alpha));
    ppll_nodes[node].thickness = thickness; // don't pack
    ppll_nodes[node].depth = depth; // don't pack
}

Node ppll_node(uint node) {
    return ppll_nodes[node];
}

uint ppll_head_node(ivec2 pixel) {
    return imageLoad(ppll_heads, pixel).r;
}

void ppll_link_node(ivec2 pixel, uint next_node) {
    uint prev_node = imageAtomicExchange(ppll_heads, pixel,
        next_node);
    ppll_nodes[next_node].prev = prev_node;
}


#endif
