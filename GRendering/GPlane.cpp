#include "GPlane.h"
#include<tuple>

GPlane::GPlane()
{
	std::tie(ID, VBO) = GMeshObject::genVABO((void *)VertexData, sizeof(VertexData));
	
	GMeshObject::EnableVertexAttrbArrays(3, 3, 2);
	TextureID =  GMeshObject::loadFromFile("Resources/block.png");
	glBindVertexArray(0);
}

void GPlane::Draw(MyShader& shader) const noexcept
{
	shader.use();
	shader.setInt("texture_diffuse", 0);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, TextureID);
	glBindVertexArray(ID);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDrawArrays(GL_TRIANGLES, 0, 6);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBindVertexArray(0);
}