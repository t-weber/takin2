/**
 * GL drawing
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 22-dec-2014
 * @license GPLv2 or GPLv3
 */

/**
 * Freetype rendering under OpenGL, inspired by:
 * https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Text_Rendering_01
 * https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Text_Rendering_02
 */

#ifndef __GL_FONT_H__
#define __GL_FONT_H__

#include "gl.h"

#define DEF_FONT "/usr/share/fonts/dejavu/DejaVuSansMono.ttf"
#define DEF_FONT_SIZE 12

#include <ft2build.h>
#include FT_FREETYPE_H

namespace tl {
class FontMap
{
	protected:
		bool m_bOk = 0;
		FT_Library m_ftLib = 0;
		FT_Face m_ftFace = 0;

		int m_iCharsPerLine=0, m_iLines=0;
		int m_iTileH=0, m_iTileW=0;
		int m_iPadH=DEF_FONT_SIZE/4, m_iPadW=DEF_FONT_SIZE/4;
		int m_iLargeW=0, m_iLargeH=0;
		unsigned char *m_pcLarge = nullptr;

		static const std::string m_strCharset;
		using t_offsmap = std::unordered_map<
			typename std::string::value_type,
			std::pair<int, int>>;
		t_offsmap m_mapOffs;

	protected:
		void UnloadFont();

	static void draw_tile(unsigned char* pcBuf,
		unsigned int iW, unsigned int iH,
		unsigned int iTileW, unsigned int iTileH,
		unsigned int iPosX, unsigned int iPosY,
		const unsigned char* pcGlyph,
		unsigned int iGlyphXOffs, unsigned int iGlyphYOffs,
		unsigned int iGlyphW, unsigned int iGlyphH);

	public:
		FontMap(const char* pcFont/*=DEF_FONT*/, int iSize/*=DEF_FONT_SIZE*/);
		FontMap();
		virtual ~FontMap();

		bool LoadFont(const char* pcFont, int iSize);
		bool LoadFont(FT_Face ftFace, int iSize);

		void dump(std::ostream& ostr) const;
		void dump(std::ostream& ostr, const std::pair<int,int>& pair) const;
		virtual bool IsOk() const { return m_bOk; }

		const unsigned char* GetBuffer() const { return m_pcLarge; }
		std::pair<int, int> GetOffset(char ch) const;

		//static std::string get_font_file(const std::string& strFind);
};


template<class T = GLdouble>
class GlFontMap : public FontMap
{
	protected:
		bool m_bOk = 0;
		GLuint m_tex = 0;
		T m_dScale = 0.01;

	protected:
		bool CreateFontTexture();

	public:
		GlFontMap() = delete;
		GlFontMap(const char* pcFont=DEF_FONT, int iSize=DEF_FONT_SIZE);
		GlFontMap(FT_Face ftFace, int iSize=DEF_FONT_SIZE);
		virtual ~GlFontMap();

		virtual bool IsOk() const override { return m_bOk && FontMap::IsOk(); }

		void BindTexture();
		void DrawText(T dX, T dY, const std::string& str, bool bCenter=1);
		void DrawText(T dX, T dY, T dZ, const std::string& str, bool bCenter=1);

		void SetScale(T dScale) { m_dScale = dScale; }
};

// --------------------------------------------------------------------------------

}

//
//#ifdef TLIBS_INC_HDR_IMPLS
//	#include "gl_font_impl.h"
//#endif

#endif
