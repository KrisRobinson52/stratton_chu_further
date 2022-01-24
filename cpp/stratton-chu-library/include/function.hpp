#pragma once

#include <array>
#include <cstddef>
#include <cstring>
#include <iostream>


/* std::function replacement to store lambda of any size on stack*/
template <size_t CtxSize, size_t Align = alignof(void*), typename Ret = void>
class function {
    alignas(Align) std::array<std::byte, CtxSize> ctx;  // lambda capture
    Ret(*func)(void*);  // operator() of lambda

    template<typename Lambda>
    void init(Ret(Lambda::*)() const, Lambda* initial_ctx)
    {
        func = [](void* ctx) -> Ret{
            auto* lambda_ptr = static_cast<Lambda*>(ctx);
            return (*lambda_ptr)();
        };
        memcpy(ctx.data(), initial_ctx, sizeof(Lambda));
    }

    template<typename Lambda>
    void init(Ret(Lambda::*)() /*mutable*/, Lambda* initial_ctx)
    {
        func = [](void* ctx) -> Ret{
            auto* lambda_ptr = static_cast<Lambda*>(ctx);
            return (*lambda_ptr)();
        };
        memcpy(ctx.data(), initial_ctx, sizeof(Lambda));
    }

public:
    template <typename Lambda>
    function& operator=(Lambda&& lmb) {
        static_assert(sizeof(Lambda) == CtxSize);
        static_assert(alignof(Lambda) == Align);

        init(&std::remove_reference_t<decltype(lmb)>::operator(), &lmb);
        return *this;
    }

    template <typename Lambda>
    function(Lambda&& lmb) {
        init(&std::remove_reference_t<decltype(lmb)>::operator(), &lmb);
    }

    function() {
        func = nullptr;
    }

    function(const function&) = default;
    function(function&&) noexcept = default;
    function& operator=(const function&) = default;
    function& operator=(function&&) = default;

    Ret operator()() {
        return func(ctx.data());
    }

    operator bool() {
        return func != nullptr;
    }
};

template <typename Lambda>
explicit function(Lambda&& lmb) ->
    function<sizeof(Lambda), alignof(Lambda), decltype(lmb())>;
