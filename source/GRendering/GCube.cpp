#include "GCube.h"

using namespace GCONST;

//TODO: �������з����ӿ�ͨ��һ��
//		�д��Ż�ʹ�� { ����ģʽ }
GCube::GCube()
{
	
	// ���ݾ������
	std::tie(model_id_, VBO) = GMeshObject::genVABO((void *)VertexData, sizeof(VertexData));

	// ���ݴ���
	GMeshObject::EnableVertexAttrbArrays(3, 3, 2);

	//TODO: ����о�������RAII�Ż�
	// Wrapper::HelperBinder����RAII
	glBindVertexArray(0);

	//TODO: ���ɴ�������Texture�ĺ�����
	// DiffuseTexureID = loadFromFile("xxx.jpg");
}

void GCube::Draw(MyShader& shader) const noexcept
{
	glBindVertexArray(model_id_);

	glDrawArrays(GL_TRIANGLES, 0, 36);
	
	glBindVertexArray(0);
}