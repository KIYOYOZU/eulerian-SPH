#include "xml_engine.h"

#include <iostream>

namespace SPH
{
//=================================================================================================//
XmlEngine::XmlEngine(const std::string &xml_name, const std::string &root_tag)
    : xml_name_(xml_name)
{
    tinyxml2::XMLElement *root = xmldoc_.NewElement(root_tag.c_str());
    xmldoc_.InsertFirstChild(root);
    root_element_ = XmlElement(root);
}
//=================================================================================================//
void XmlEngine::addElementToXmlDoc(const std::string &element_name)
{
    tinyxml2::XMLElement *el = xmldoc_.NewElement(element_name.c_str());
    root_element_.ptr_->InsertEndChild(el);
}
//=================================================================================================//
void XmlEngine::addChildToElement(XmlElement &father_element, const std::string &child_name)
{
    if (!father_element.ptr_)
        return;
    tinyxml2::XMLElement *child = xmldoc_.NewElement(child_name.c_str());
    father_element.ptr_->InsertEndChild(child);
}
//=================================================================================================//
XmlElement XmlEngine::getChildElement(const std::string &tag)
{
    if (!root_element_.ptr_)
        return XmlElement(nullptr);
    return XmlElement(root_element_.ptr_->FirstChildElement(tag.c_str()));
}
//=================================================================================================//
void XmlEngine::writeToXmlFile(const std::string &filefullpath)
{
    tinyxml2::XMLError result = xmldoc_.SaveFile(filefullpath.c_str());
    if (result != tinyxml2::XML_SUCCESS)
    {
        std::cerr << "XmlEngine::writeToXmlFile failed for: " << filefullpath
                  << "  error=" << result << std::endl;
        std::exit(1);
    }
}
//=================================================================================================//
void XmlEngine::loadXmlFile(const std::string &filefullpath)
{
    tinyxml2::XMLError result = xmldoc_.LoadFile(filefullpath.c_str());
    if (result != tinyxml2::XML_SUCCESS)
    {
        std::cerr << "XmlEngine::loadXmlFile failed for: " << filefullpath
                  << "  error=" << result << std::endl;
        std::exit(1);
    }
    root_element_ = XmlElement(xmldoc_.RootElement());
}
//=================================================================================================//
std::string XmlEngine::getRootElementTag()
{
    if (!root_element_.ptr_)
        return "";
    const char *name = root_element_.ptr_->Name();
    return name ? std::string(name) : "";
}
//=================================================================================================//
std::string XmlEngine::getElementTag(XmlElement &element)
{
    return element.getElementTag();
}
//=================================================================================================//
void XmlEngine::resizeXmlDocForParticles(size_t input_size)
{
    size_t total_elements = root_element_.childElementCount();
    if (total_elements <= input_size)
    {
        for (size_t i = total_elements; i != input_size; ++i)
            addElementToXmlDoc("particle");
    }
    else
    {
        std::cout << "\n Error: XML Engine allows increase data size only!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        std::exit(1);
    }
}
//=================================================================================================//
size_t XmlEngine::SizeOfXmlDoc()
{
    return root_element_.childElementCount();
}
//=================================================================================================//
} // namespace SPH
