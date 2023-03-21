#pragma once
#include "MyShader.h"
#include <stb_image.h>

namespace GCONST {
	constexpr size_t FLOAT_SIZE = sizeof(float);
	constexpr size_t INT_SIZE = sizeof(int);

}

namespace WARRPER {
	struct HelperBinder{
		HelperBinder() {};
		~HelperBinder() { glEnableVertexAttribArray(0); }
	};
}

namespace GComponent{
	using namespace glm;
	using std::vector;
	struct Vertex {
		vec3 pos;
		vec3 norm;
		vec2 coord;
	};
	struct Triangle {
		int first;
		int second;
		int third;
	};

	struct Line {
		int first;
		int second;
	};

	class MeshComponent {
	public:
		vector<Vertex> vertices;
		vector<Triangle> indicies;
		MeshComponent()  = default;
		~MeshComponent() = default;
		MeshComponent(MeshComponent& m) = default;
		MeshComponent(MeshComponent&& m) = default;
		MeshComponent& operator=(MeshComponent& m) = default;
		MeshComponent& operator=(MeshComponent&& m) = default;
		MeshComponent(vector<Vertex> vertices, vector<unsigned> indicies);
	};
}

class GMeshObject
{
public:
	unsigned model_id_;
	unsigned VBO;

	GMeshObject() = default;

	virtual ~GMeshObject() = 0 { glDeleteBuffers(1, &VBO); glDeleteVertexArrays(1, &model_id_); };

	virtual void Draw(MyShader& shader) const noexcept= 0 {};

#pragma region StaticFunc

	/// <summary>
	/// ���ļ���ȡͼƬ������unsigned��ŵ�����
	/// </summary>
	/// <param name="fileName"> {const string&} �ļ��������ӣ�image.png </param>
	/// <param name="repeat"> {bool} �ظ���־��ͼƬ��ȡ���ȡ�Ƿ���� </param>
	/// <returns> {unsigned} ���ؾ����ָ���ڲ�������洢λ�� </returns>
	static unsigned 
	loadFromFile(const std::string& fileName, bool repeat = true)
	// Ŀǰ������ ID δʹ������ڶ�ʧ����
	{
		unsigned tempID;
		glGenTextures(1, &tempID);
		int width, height, nrChannels;
		unsigned char* data = stbi_load(fileName.c_str(), &width, &height, &nrChannels, 0);
		if (data != nullptr)
		{
			GLenum format = 0;
			switch (nrChannels)
			{
				case 1:format = GL_RED; break;
				case 3:format = GL_RGB; break;
				case 4:format = GL_RGBA;
			}
			glBindTexture(GL_TEXTURE_2D, tempID);
			glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
			glGenerateMipmap(GL_TEXTURE_2D);
			GLint wrapMethod;
			if (repeat)
			{
				wrapMethod = GL_REPEAT;
			}
			else
			{
				wrapMethod = GL_CLAMP_TO_EDGE;
			}
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapMethod);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapMethod);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

			stbi_image_free(data);
		}
		else
		{
			std::cout << "GObject::loadFromFile::ERROR::name->" << fileName << std::endl;
		}
		return tempID;
	}

	/// <summary>
	/// ���ɾ���ӿڿ黺�棬�����ھ������Ⲽ�ֵ����
	/// </summary>
	/// <returns> {unsigned} ���ؾ����ָ��ӿڿ���ڲ�����λ��</returns>
	static unsigned 
	genMatrices() {
		unsigned uboMatrices;
		glGenBuffers(1, &uboMatrices);
		glBindBuffer(GL_UNIFORM_BUFFER, uboMatrices);
		glBufferData(GL_UNIFORM_BUFFER, 2 * sizeof(glm::mat4), NULL, GL_STATIC_DRAW);
		glBindBufferBase(GL_UNIFORM_BUFFER, 0, uboMatrices);
		return uboMatrices;
	}
	
	/// <summary>
	/// ���þ����ӿڵĲ��������������Ⲽ�֣����ڴ��ݹ��õĶ����������ξ���
	/// </summary>
	/// <param name="ubo">  {unsigned}  buffer��ID��ָ��ӿڿ��λ��								</param>
	/// <param name="pro">  {glm::mat4} projection ������������ѹ���ӿ�							</param>
	/// <param name="view"> {glm::mat4} view	   �ӿھ������ڰ��������������ϵת�����������ϵ </param>
	static void 
	setMatrices(unsigned ubo, glm::mat4 pro, glm::mat4 view)
	{
		glBindBuffer(GL_UNIFORM_BUFFER, ubo);
		glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(glm::mat4), glm::value_ptr(pro));
		glBufferSubData(GL_UNIFORM_BUFFER, sizeof(glm::mat4), sizeof(glm::mat4), glm::value_ptr(view));
		glBindBuffer(GL_UNIFORM_BUFFER, 0);
	}

	/// <summary>
	/// �����ݰ��������ڴ���ϣ��������飬���ڰ󶨲�����VertexArray
	/// </summary>
	/// <param name="data">  {float[]} ���飬���ڴ����贫��Shader������				</param>
	/// <param name="size">  {int}	   ��С�����ڴ������ݵĴ�С����ֹԽ��			</param>
	/// <returns> {tuple<u, u>} ���ؾ���� ���ڷ��ض���λ���Լ��������ݻ����λ��	</returns>
	static std::pair<unsigned, unsigned>
		genVABO(void* data, size_t size)
	{
		unsigned VAO, VBO;
		glGenVertexArrays(1, &VAO);
		glGenBuffers(1, &VBO);
		glBindVertexArray(VAO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, size, data, GL_STATIC_DRAW);
		return { VAO, VBO };
	}

	/// <summary>
	/// ������������Ű���˳��������;���Լ����巵��BufferIdx
	/// </summary>
	/// <param name="data">  {float[]} ���飬���ڴ����贫��Shader������				</param>
	/// <param name="size">  {int}	   ��С�����ڴ������ݵĴ�С����ֹԽ��			</param>
	/// <returns> {unsigned} ���ؾ���� ���ڷ���Ԫ������˳��						</returns>
	template<typename... Args>
	static 
	unsigned genEBO(Args & ... _args)
	{
		static_assert(std::conjunction_v<std::_Is_specialization<Args, std::vector>...>,
			"The Input Type can't cast Vector, ERROR::GObject::genEBO");
		unsigned EBO;
		glGenBuffers(1, &EBO);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		size_t size = ((sizeof(std::remove_reference_t<Args>::value_type) * _args.size()) +...);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, nullptr, GL_STATIC_DRAW);
		{
			size_t offset = 0, tempSize = 0;
			((tempSize = sizeof(std::remove_reference_t<Args>::value_type) * _args.size(),
			  glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, offset, tempSize, &_args[0]),
			  offset += tempSize),
			  ...);
		}
		return EBO;
	}

	/// <summary>
	/// ���ڼ��������ָ�룬 ʹ�����ܹ����ݵ���ȷ��ָ��λ��
	/// ��ע�⡿�� ���ִ����һ������ͣ� �������ʹ�õ������
	/// </summary>
	/// <typeparam name="...Args">������</typeparam>
	/// <param name="..._args"></param>
	template<class ... Args>
	static void EnableVertexAttrbArrays(Args... _args)
	{
		static_assert(std::conjunction_v<std::is_same<int, Args>...>,
			"The Input Type is not Int, ERROR::GObject::EnableVertexAttribArrays");
		const unsigned TOTAL = (_args + ...);
		unsigned idx = 0, loc = 0;
		
		((glEnableVertexAttribArray(idx), 
		  glVertexAttribPointer(idx++, 
							  _args, 
			                  GL_FLOAT, 
			                  GL_FALSE, 
						      TOTAL * GCONST::FLOAT_SIZE,
			                  (void*)(loc * GCONST::FLOAT_SIZE)), 
			loc += _args), ...);
	}
	
	template<class... Args>
	static void EnableVertexAttribArrays_continus(Args&... _args)
	{
		static_assert(std::conjunction_v<std::_Is_specialization<Args, std::vector>...>,
			"The Input Type can't cast Vector, ERROR::GObject::genEBO");
		unsigned idx = 0, stride = 0;
		size_t singleSize = 0, loc = 0;
		
		((  singleSize = sizeof(std::remove_reference_t<Args>::value_type),
			stride =  singleSize/ sizeof(float),
			glEnableVertexAttribArray(idx),
			glVertexAttribPointer(idx++, stride, GL_FLOAT, GL_FALSE, singleSize, (void*)(loc)),
			loc += singleSize * _args.size()), ...);

	}
#pragma endregion
	
};

