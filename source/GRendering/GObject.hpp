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
	/// 从文件读取图片并生成unsigned标号的纹理
	/// </summary>
	/// <param name="fileName"> {const string&} 文件名，例子：image.png </param>
	/// <param name="repeat"> {bool} 重复标志，图片读取后采取是否包覆 </param>
	/// <returns> {unsigned} 返回句柄，指向内部该纹理存储位置 </returns>
	static unsigned 
	loadFromFile(const std::string& fileName, bool repeat = true)
	// 目前若返回 ID 未使用则存在丢失风险
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
	/// 生成矩阵接口块缓存，适用于具有特殊布局的情况
	/// </summary>
	/// <returns> {unsigned} 返回句柄，指向接口块的内部输入位置</returns>
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
	/// 设置矩阵块接口的参数，适用于特殊布局，用于传递公用的二个相机内外参矩阵
	/// </summary>
	/// <param name="ubo">  {unsigned}  buffer的ID，指向接口块的位置								</param>
	/// <param name="pro">  {glm::mat4} projection 放缩矩阵，用于压缩视口							</param>
	/// <param name="view"> {glm::mat4} view	   视口矩阵，用于把物体从世界坐标系转换到相机坐标系 </param>
	static void 
	setMatrices(unsigned ubo, glm::mat4 pro, glm::mat4 view)
	{
		glBindBuffer(GL_UNIFORM_BUFFER, ubo);
		glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(glm::mat4), glm::value_ptr(pro));
		glBufferSubData(GL_UNIFORM_BUFFER, sizeof(glm::mat4), sizeof(glm::mat4), glm::value_ptr(view));
		glBindBuffer(GL_UNIFORM_BUFFER, 0);
	}

	/// <summary>
	/// 将数据绑定在两个内存块上，传入数组，用于绑定并返回VertexArray
	/// </summary>
	/// <param name="data">  {float[]} 数组，用于传递需传给Shader的数组				</param>
	/// <param name="size">  {int}	   大小，用于传递数据的大小，防止越界			</param>
	/// <returns> {tuple<u, u>} 返回句柄， 用于返回顶点位置以及顶点数据缓存的位置	</returns>
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
	/// 将顶点连接序号按，顺序打出，用途由自己定义返回BufferIdx
	/// </summary>
	/// <param name="data">  {float[]} 数组，用于传递需传给Shader的数组				</param>
	/// <param name="size">  {int}	   大小，用于传递数据的大小，防止越界			</param>
	/// <returns> {unsigned} 返回句柄， 用于返回元素连接顺序						</returns>
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
	/// 用于激活传递数据指针， 使数据能够传递到正确的指针位置
	/// 【注意】： 比手打的慢一次总求和， 但提高了使用的灵活性
	/// </summary>
	/// <typeparam name="...Args">参数包</typeparam>
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

