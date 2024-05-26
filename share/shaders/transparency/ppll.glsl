#ifndef VKHR_PPLL_GLSL
#define VKHR_PPLL_GLSL

#define PPLL_NULL_NODE 0xffffffff
#define PPLL_MAXIMUM_DEPTH 3.e+38

struct Node {
    uint color;
    float depth;
    uint prev;
};

layout(binding = 5, r32ui) uniform uimage2D ppll_heads;
layout(binding = 6, std430) buffer LinkedList {
    Node ppll_nodes[];
};

layout(binding = 7) uniform Config { uint ppll_size; };
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

void ppll_node_data(uint node, vec4 color, float depth) {
    ppll_nodes[node].color = packUnorm4x8(color);
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

void ppll_link_node_sorted2(ivec2 pixel, uint node) {
    uint head = ppll_head_node(pixel);
    float depth = ppll_nodes[node].depth;

    if (ppll_nodes[head].depth < depth) {
        ppll_link_node(pixel, node);
        return;
    }

    uint nextNode = head;
    uint prevNode = ppll_nodes[head].prev;

    while (prevNode != PPLL_NULL_NODE && ppll_nodes[prevNode].depth > depth) {
        nextNode = prevNode;
        prevNode = ppll_nodes[prevNode].prev;

    }

    ppll_nodes[nextNode].prev = node;
    ppll_nodes[node].prev = prevNode;
}

bool fragment_relevant(ivec2 pixel, float depth, float eps) {
    uint node_ind = ppll_head_node(pixel);
    while (node_ind != PPLL_NULL_NODE) {
        Node node = ppll_node(node_ind);
        float node_depth = node.depth;
        if (node_depth - eps < depth && depth < node_depth + eps) return false;
        node_ind = node.prev;
    }
    return true;
}

bool fragment_relevant_sorted(ivec2 pixel, float depth, float eps) {
    uint node_ind = ppll_head_node(pixel);
    while (node_ind != PPLL_NULL_NODE) {
        Node node = ppll_node(node_ind);
        float node_depth = node.depth;
        if (node_depth - eps < depth && depth < node_depth + eps) return false;
        if (node_depth - eps > depth) return true;
        node_ind = node.prev;
    }
    return true;
}
#endif
