
// Derived from: https://en.wikipedia.org/wiki/MurmurHash

static inline uint32_t murmur_32_scramble(uint32_t k) {
    k *= 0xcc9e2d51;
    k = (k << 15) | (k >> 17);
    k *= 0x1b873593;
    return k;
}

uint32_t murmur3_vec(uint32_t *data, uint32_t size, uint32_t seed) {
    uint32_t h = seed;
    for (uint32_t i = 0; i < size; i++) {
        h ^= murmur_32_scramble(data[i]);
        h = (h << 13) | (h >> 19);
        h = h * 5 + 0xe6546b64;
    }
    h ^= size;
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}
