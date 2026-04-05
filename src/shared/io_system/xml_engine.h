/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * xml_engine.h - XML engine based on tinyxml2 (no Simbody dependency).     *
 * Provides XmlEngine class with interface compatible with the original      *
 * SimTK-based implementation.                                               *
 * ------------------------------------------------------------------------- */
/**
 * @file    xml_engine.h
 * @brief   XML engine built on tinyxml2 for xml input/output.
 * @author  Adapted from Chi Zhang and Xiangyu Hu (original SimTK version)
 */

#ifndef XML_ENGINE_H
#define XML_ENGINE_H

#include "base_data_type_package.h"
#include "tinyxml2.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;

namespace SPH
{

/**
 * @class XmlElement
 * @brief Wrapper around tinyxml2::XMLElement* that mimics SimTK::Xml::Element interface.
 *
 * This allows existing code using SimTK::Xml::Element variables to work
 * without modification by re-mapping the type alias below.
 */
class XmlElement
{
  public:
    tinyxml2::XMLElement *ptr_ = nullptr;

    XmlElement() = default;
    explicit XmlElement(tinyxml2::XMLElement *p) : ptr_(p) {}

    /** Returns true if the underlying pointer is valid. */
    bool isValid() const { return ptr_ != nullptr; }

    /**
     * Returns an iterator (pointer) to the first child element,
     * optionally matching a tag name.
     */
    tinyxml2::XMLElement *element_begin(const std::string &tag = "") const
    {
        if (!ptr_)
            return nullptr;
        if (tag.empty())
            return ptr_->FirstChildElement();
        return ptr_->FirstChildElement(tag.c_str());
    }

    /** Returns the sentinel end iterator (nullptr). */
    tinyxml2::XMLElement *element_end() const
    {
        return nullptr;
    }

    /** Appends a new child element with the given name. */
    void insertChildElement(const std::string &child_name)
    {
        if (!ptr_)
            return;
        tinyxml2::XMLDocument *doc = ptr_->GetDocument();
        tinyxml2::XMLElement *child = doc->NewElement(child_name.c_str());
        ptr_->InsertEndChild(child);
    }

    /** Returns the tag name of this element. */
    std::string getElementTag() const
    {
        if (!ptr_)
            return "";
        const char *name = ptr_->Name();
        return name ? std::string(name) : "";
    }

    /** Returns an XmlElement for the first child matching tag (may be invalid). */
    XmlElement getOptionalElement(const std::string &tag) const
    {
        if (!ptr_)
            return XmlElement(nullptr);
        return XmlElement(ptr_->FirstChildElement(tag.c_str()));
    }

    /** Compute the number of direct child elements. */
    std::size_t childElementCount() const
    {
        std::size_t count = 0;
        for (tinyxml2::XMLElement *e = element_begin(); e != element_end(); e = e->NextSiblingElement())
            ++count;
        return count;
    }
};

/**
 * @class XmlEngine
 * @brief XML document engine backed by tinyxml2.
 *
 * Public interface is designed to be a drop-in replacement for the original
 * SimTK-based XmlEngine so that call-sites require minimal changes.
 */
class XmlEngine
{
  protected:
    std::string xml_name_;
    tinyxml2::XMLDocument xmldoc_;

  public:
    /** Public root element (wraps the tinyxml2 root element pointer). */
    XmlElement root_element_;

    /** Iterator type alias — kept compatible with old element_iterator usage. */
    using element_iterator = tinyxml2::XMLElement *;

    XmlEngine(const std::string &xml_name, const std::string &root_tag);
    virtual ~XmlEngine() {}

    /** Add a new child element directly under the root. */
    void addElementToXmlDoc(const std::string &element_name);

    /** Add a child element to an arbitrary parent XmlElement. */
    void addChildToElement(XmlElement &father_element, const std::string &child_name);

    /** Get the first direct child of root with the given tag. */
    XmlElement getChildElement(const std::string &tag);

    /** Write the XML document to file. */
    void writeToXmlFile(const std::string &filefullpath);

    /** Load/parse an XML file. */
    void loadXmlFile(const std::string &filefullpath);

    /** Get root tag string. */
    std::string getRootElementTag();

    /** Get element tag string. */
    std::string getElementTag(XmlElement &element);

    /** Resize (by adding) child elements to reach input_size. */
    void resizeXmlDocForParticles(size_t input_size);

    /** Count of direct children of root. */
    size_t SizeOfXmlDoc();

    //----------------------------------------------------------------------
    // setAttributeToElement — scalar types
    //----------------------------------------------------------------------
    template <typename T>
    void setAttributeToElement(element_iterator ele_ite,
                               const std::string &attrib_name,
                               const T &value)
    {
        if (!ele_ite)
            return;
        std::ostringstream oss;
        oss << value;
        ele_ite->SetAttribute(attrib_name.c_str(), oss.str().c_str());
    }

    //----------------------------------------------------------------------
    // setAttributeToElement — Eigen column vector
    //----------------------------------------------------------------------
    template <int DIMENSION, auto... Rest>
    void setAttributeToElement(element_iterator ele_ite,
                               const std::string &attrib_name,
                               const Eigen::Matrix<Real, DIMENSION, 1, Rest...> &value)
    {
        if (!ele_ite)
            return;
        std::ostringstream oss;
        oss << "~[";
        for (int i = 0; i < DIMENSION; ++i)
        {
            if (i)
                oss << ",";
            oss << value[i];
        }
        oss << "]";
        ele_ite->SetAttribute(attrib_name.c_str(), oss.str().c_str());
    }

    //----------------------------------------------------------------------
    // setAttributeToElement — Eigen square matrix
    //----------------------------------------------------------------------
    template <int DIMENSION, auto... Rest>
    void setAttributeToElement(element_iterator ele_ite,
                               const std::string &attrib_name,
                               const Eigen::Matrix<Real, DIMENSION, DIMENSION, Rest...> &value)
    {
        if (!ele_ite)
            return;
        std::ostringstream oss;
        oss << "~[";
        for (int r = 0; r < DIMENSION; ++r)
        {
            if (r)
                oss << ";";
            oss << "(";
            for (int c = 0; c < DIMENSION; ++c)
            {
                if (c)
                    oss << ",";
                oss << value(r, c);
            }
            oss << ")";
        }
        oss << "]";
        ele_ite->SetAttribute(attrib_name.c_str(), oss.str().c_str());
    }

    //----------------------------------------------------------------------
    // getRequiredAttributeValue — scalar types
    //----------------------------------------------------------------------
    template <typename T>
    void getRequiredAttributeValue(element_iterator ele_ite_,
                                   const std::string &attrib_name,
                                   T &value)
    {
        if (!ele_ite_)
            return;
        const char *attr = ele_ite_->Attribute(attrib_name.c_str());
        if (!attr)
        {
            std::cerr << "XmlEngine: attribute '" << attrib_name << "' not found." << std::endl;
            std::exit(1);
        }
        std::istringstream iss(attr);
        iss >> value;
    }

    //----------------------------------------------------------------------
    // getRequiredAttributeValue — Eigen column vector
    //----------------------------------------------------------------------
    template <int DIMENSION, auto... Rest>
    void getRequiredAttributeValue(element_iterator ele_ite_,
                                   const std::string &attrib_name,
                                   Eigen::Matrix<Real, DIMENSION, 1, Rest...> &value)
    {
        if (!ele_ite_)
            return;
        const char *attr = ele_ite_->Attribute(attrib_name.c_str());
        if (!attr)
        {
            std::cerr << "XmlEngine: attribute '" << attrib_name << "' not found." << std::endl;
            std::exit(1);
        }
        // Format written by setAttributeToElement: ~[v0,v1,...,vN-1]
        std::string s(attr);
        // Strip leading "~[" and trailing "]"
        auto start = s.find('[');
        auto end = s.rfind(']');
        if (start == std::string::npos || end == std::string::npos)
        {
            std::cerr << "XmlEngine: malformed vector attribute '" << attrib_name << "': " << s << std::endl;
            std::exit(1);
        }
        std::string inner = s.substr(start + 1, end - start - 1);
        std::istringstream iss(inner);
        std::string token;
        for (int i = 0; i < DIMENSION; ++i)
        {
            if (!std::getline(iss, token, ','))
                break;
            value[i] = static_cast<Real>(std::stod(token));
        }
    }

    //----------------------------------------------------------------------
    // getRequiredAttributeValue — Eigen square matrix
    //----------------------------------------------------------------------
    template <int DIMENSION, auto... Rest>
    void getRequiredAttributeValue(element_iterator ele_ite_,
                                   const std::string &attrib_name,
                                   Eigen::Matrix<Real, DIMENSION, DIMENSION, Rest...> &value)
    {
        if (!ele_ite_)
            return;
        const char *attr = ele_ite_->Attribute(attrib_name.c_str());
        if (!attr)
        {
            std::cerr << "XmlEngine: attribute '" << attrib_name << "' not found." << std::endl;
            std::exit(1);
        }
        // Format: ~[(r0c0,r0c1,...);(r1c0,...);...]
        std::string s(attr);
        auto start = s.find('[');
        auto end = s.rfind(']');
        if (start == std::string::npos || end == std::string::npos)
        {
            std::cerr << "XmlEngine: malformed matrix attribute '" << attrib_name << "': " << s << std::endl;
            std::exit(1);
        }
        std::string inner = s.substr(start + 1, end - start - 1);
        // Split rows by ';'
        std::istringstream row_stream(inner);
        std::string row_str;
        for (int r = 0; r < DIMENSION; ++r)
        {
            if (!std::getline(row_stream, row_str, ';'))
                break;
            // Remove surrounding parentheses
            auto ps = row_str.find('(');
            auto pe = row_str.rfind(')');
            std::string row_inner = (ps != std::string::npos && pe != std::string::npos)
                                        ? row_str.substr(ps + 1, pe - ps - 1)
                                        : row_str;
            std::istringstream col_stream(row_inner);
            std::string token;
            for (int c = 0; c < DIMENSION; ++c)
            {
                if (!std::getline(col_stream, token, ','))
                    break;
                value(r, c) = static_cast<Real>(std::stod(token));
            }
        }
    }
};

} // namespace SPH

// ---------------------------------------------------------------------------
// Compatibility shim: map SimTK::Xml::Element and SimTK::Xml::element_iterator
// to our tinyxml2-based replacements so that call-site code requires no change.
// ---------------------------------------------------------------------------
namespace SimTK
{

/**
 * @brief Minimal SpatialVec shim — replaces SimTK::SpatialVec (rotation + force 6-vector).
 * Layout: [0] = angular (torque/rotation), [1] = linear (force/translation).
 */
struct SpatialVec
{
    Eigen::Matrix<double, 3, 1> data_[2];

    SpatialVec() { data_[0].setZero(); data_[1].setZero(); }
    SpatialVec(const Eigen::Matrix<double, 3, 1> &ang,
               const Eigen::Matrix<double, 3, 1> &lin)
    { data_[0] = ang; data_[1] = lin; }

    Eigen::Matrix<double, 3, 1> &operator[](int i) { return data_[i]; }
    const Eigen::Matrix<double, 3, 1> &operator[](int i) const { return data_[i]; }
};

namespace Xml
{
using Element = SPH::XmlElement;
using element_iterator = tinyxml2::XMLElement *;
} // namespace Xml
} // namespace SimTK

#endif // XML_ENGINE_H
