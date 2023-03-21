#pragma once
#include "GObject.hpp"

class GPlane :
    public GMeshObject
{
public:
    GPlane();
    ~GPlane() { glDeleteVertexArrays(1, &model_id_); glDeleteBuffers(1, &VBO);}
    void Draw(MyShader& shader, unsigned TexID) const {};
    void Draw(MyShader& shader) const noexcept override;
private:
    static constexpr float VertexData[] = {
        // positions                            // texture Coords (note we set these higher than 1 (together with GL_REPEAT as texture wrapping mode). this will cause the floor texture to repeat)
         5.0f, -0.5f,  5.0f,  0.0f, 1.0f, 0.0f,  8.0f, 0.0f,
        -5.0f, -0.5f, -5.0f,  0.0f, 1.0f, 0.0f,  0.0f, 8.0f,
        -5.0f, -0.5f,  5.0f,  0.0f, 1.0f, 0.0f,  0.0f, 0.0f,

         5.0f, -0.5f,  5.0f,  0.0f, 1.0f, 0.0f,  8.0f, 0.0f,
         5.0f, -0.5f, -5.0f,  0.0f, 1.0f, 0.0f,  8.0f, 8.0f,
        -5.0f, -0.5f, -5.0f,  0.0f, 1.0f, 0.0f,  0.0f, 8.0f,
    };
    unsigned TextureID;
};

