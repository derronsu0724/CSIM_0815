/*
 *  FastCSIM Copyright (C) 2021 cassuto                                    
 *  This project is free edition; you can redistribute it and/or           
 *  modify it under the terms of the GNU Lesser General Public             
 *  License(LGPL) as published by the Free Software Foundation; either      
 *  version 2.1 of the License, or (at your option) any later version.     
 *                                                                         
 *  This project is distributed in the hope that it will be useful,        
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU      
 *  Lesser General Public License for more details.                        
 */

#ifndef CSIM_MODELLOADER_H_
#define CSIM_MODELLOADER_H_

#include <map>
#include "csim/model/ModelBase.h"

// 定义命名空间csimModel
namespace csimModel
{
    class ModelBase;
}

// 定义命名空间csim
namespace csim
{
    class Circuit;

    // 定义模型条目类
    class ModelEntry
    {
    public:
        // 构造函数
        ModelEntry(void *handle, csimModel::pfnCreateModel_t pfnCreate, csimModel::pfnDeleteModel_t pfnDelete,
                   const ModelDescriptor *descriptor,
                   const PropertyMdlDescriptor *mdlDescriptors,
                   size_t numMdlDescriptors);
        // 析构函数
        ~ModelEntry();

        // 创建模型实例
        csimModel::ModelBase *createInstance(MODELBASE_CONSTRUCTOR_DEF) const;

        // 删除模型实例
        void deleteInstance(csimModel::ModelBase *model) const;

        // 获取模型描述符
        inline const ModelDescriptor *descriptor() const
        {
            return m_descriptor;
        }

        // 获取模型描述符数量
        inline size_t numMdlDescriptors() const
        {
            return m_numMdlDescriptors;
        }

        // 获取指定索引的模型描述符
        inline const PropertyMdlDescriptor *mdlDescriptors(size_t i) const
        {
            return m_mdlDescriptors + i;
        }

        // 定义模型条目类
        class MdlEntry
        {
        public:
            // 创建模型实例
            csimModel::PropertyMdl *createInstance() const;
            // 删除模型实例
            void deleteInstance(csimModel::PropertyMdl *mdl) const;
            // 获取属性描述符
            const PropertyMdlPropDescriptor *getProperty(const char *name) const;

        public:
            // 模型描述符
            const PropertyMdlDescriptor *desc;
            // 属性描述符映射
            std::map<std::string, const PropertyMdlPropDescriptor *> props;
        };

        // 获取模型条目
        int getMdlEntry(const char *mdl, const MdlEntry **out) const;

    private:
        // 构建模型索引
        int buildMdlIndexs();

    private:
        // 动态链接库句柄
        void *m_dllHandle;
        // 创建模型函数指针
        csimModel::pfnCreateModel_t m_pfnCreate;
        // 删除模型函数指针
        csimModel::pfnDeleteModel_t m_pfnDelete;
        // 模型描述符
        const ModelDescriptor *m_descriptor;
        // 模型描述符数组
        const PropertyMdlDescriptor *m_mdlDescriptors;
        // 模型描述符数量
        size_t m_numMdlDescriptors;
        // 模型映射
        std::map<std::string, MdlEntry> m_mdls;
    };

    // 定义模型加载器类
    class ModelLoader
    {
    public:
        // 加载模型
        static ModelEntry *load(const char *filename);
    };

}

#endif // CSIM_MODELLOADER_H_