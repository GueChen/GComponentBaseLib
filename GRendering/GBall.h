#pragma once
#include "GObject.hpp"

namespace PhsycalLagrangian {
    using namespace glm;
    using namespace GComponent;
    class GBall :
        public GMeshObject
    {
    public:
        GBall() = delete;
        GBall(vec3 o, float r);
        ~GBall() {}
        void setupMesh(int lo, int la);
        void Draw(MyShader& shader) const noexcept override;
        bool CalcCollision(vec3& p) const { return length(p - center) - 0.01f < radius; }
        vec3 GetOutDirection(vec3& p) const { return  center + (radius + 0.01f) * (p - center)/length(p - center); }
        void MoveBall(const vec3& v) { center += v; }
    private:
        int EBO;
        unsigned TextureID;
        vec3 center;
        float radius;
        Mesh mesh;
    };
}

