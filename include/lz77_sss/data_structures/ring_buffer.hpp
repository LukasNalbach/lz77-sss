#pragma once

#include <bit>

template <typename val_t, typename pos_t>
class ring_buffer {
    protected:
    std::vector<val_t> data;
    const pos_t buffer_size;
    const pos_t mod_mask;
    pos_t virtual_size;
    
    public:
    ring_buffer(pos_t buffer_size) :
        buffer_size(std::bit_ceil(buffer_size)),
        mod_mask(buffer_size - 1),
        virtual_size(0),
        data(buffer_size) {}

    inline void emplace_back(val_t&& val) {
        data[virtual_size++ & mod_mask] = std::move(val);
    }

    inline pos_t size() const {
        return virtual_size;
    }

    inline void resize(pos_t size) {
        virtual_size = size;
    }

    inline void advance() {
        virtual_size++;
    }

    inline val_t& operator[](pos_t const index) {
        return data[index & mod_mask];
    }

    inline val_t& front() {
        return operator[](0);
    }

    inline val_t& back() {
        return operator[](virtual_size - 1);
    }
};