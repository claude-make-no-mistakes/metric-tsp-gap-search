#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

// Calls originate from search.py and verify.py. Vertex 0 anchors the tour
// and is excluded from the enumerated side of each cut.
extern "C" int64_t minimum_cut(int n, const int64_t* edges, int* side) {
    int64_t degree[20] = {};
    int64_t incident[20] = {};
    for (int u = 0; u < n; ++u) {
        for (int v = 0; v < n; ++v) {
            degree[u] += edges[u * n + v];
        }
    }
    int64_t cut = 0;
    int64_t best = std::numeric_limits<int64_t>::max();
    const int end = 1 << (n - 1);
    int previous = 0;
    for (int index = 1; index < end; ++index) {
        const int mask = index ^ (index >> 1);
        const int vertex = __builtin_ctz(mask ^ previous) + 1;
        const int sign = (mask > previous) ? 1 : -1;
        cut += sign * (degree[vertex] - 2 * incident[vertex]);
        for (int v = 1; v < n; ++v) {
            incident[v] += sign * edges[vertex * n + v];
        }
        if (cut < best) {
            best = cut;
            *side = mask << 1;
        }
        previous = mask;
    }
    return best;
}

extern "C" int64_t shortest_tour(int n, const int64_t* distance, int* tour) {
    const int width = n - 1;
    const int end = 1 << width;
    const int64_t infinity = std::numeric_limits<int64_t>::max() / 4;
    std::vector<int64_t> path(end * width, infinity);
    for (int v = 0; v < width; ++v) {
        path[(1 << v) * width + v] = distance[v + 1];
    }
    for (int mask = 1; mask < end; ++mask) {
        for (int vertices = mask; vertices; vertices &= vertices - 1) {
            const int v = __builtin_ctz(vertices);
            const int remaining = mask ^ (1 << v);
            const int64_t* previous = &path[remaining * width];
            auto& cost = path[mask * width + v];
            for (int choices = remaining; choices; choices &= choices - 1) {
                const int u = __builtin_ctz(choices);
                cost = std::min(cost, previous[u] + distance[(u + 1) * n + v + 1]);
            }
        }
    }
    int mask = end - 1;
    int endpoint = 0;
    int64_t best = infinity;
    for (int v = 0; v < width; ++v) {
        const int64_t cost = path[mask * width + v] + distance[(v + 1) * n];
        if (cost < best) {
            best = cost;
            endpoint = v;
        }
    }
    tour[0] = 0;
    for (int position = width; position > 0; --position) {
        tour[position] = endpoint + 1;
        const int64_t cost = path[mask * width + endpoint];
        mask ^= 1 << endpoint;
        for (int choices = mask; choices; choices &= choices - 1) {
            const int u = __builtin_ctz(choices);
            if (path[mask * width + u] + distance[(u + 1) * n + endpoint + 1] == cost) {
                endpoint = u;
                break;
            }
        }
    }
    return best;
}
