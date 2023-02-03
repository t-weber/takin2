/**
 * GL drawing
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 22-dec-2014
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

/**
 * Freetype rendering under OpenGL, inspired by:
 * https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Text_Rendering_01
 * https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Text_Rendering_02
 */

#ifndef __GL_FONT_HELPER_IMPL_H__
#define __GL_FONT_HELPER_IMPL_H__

#include "gl_font.h"
#include "../log/log.h"
#include "../helper/misc.h"
#include "../string/string.h"

#include <glext.h>
#include <iomanip>


// --------------------------------------------------------------------------------

namespace tl {

void FontMap::draw_tile(unsigned char* pcBuf,
	unsigned int iW, unsigned int iH,
	unsigned int iTileW, unsigned int iTileH,
	unsigned int iPosX, unsigned int iPosY,
	const unsigned char* pcGlyph,
	unsigned int iGlyphXOffs, unsigned int iGlyphYOffs,
	unsigned int iGlyphW, unsigned int iGlyphH)
{
	unsigned int iYOffs = iPosY*iTileH + iGlyphYOffs;
	unsigned int iXOffs = iPosX*iTileW + iGlyphXOffs;

	for(unsigned int iY=0; iY<iGlyphH; ++iY)
	{
		unsigned int iYLarge = iY+iYOffs;
		if(iYLarge >= iH) break;

		if(iYLarge >= (iPosY+1)*iTileH)
		{
			break;
		}

		for(unsigned int iX=0; iX<iGlyphW; ++iX)
		{
			unsigned int iXLarge = iX+iXOffs;
			if(iXLarge >= iW) break;

			if(iXLarge >= (iPosX+1)*iTileW)
			{
				break;
			}


			pcBuf[iYLarge*iW + iXLarge] = pcGlyph[iY*iGlyphW + iX];

		}
	}
}

const std::string FontMap::m_strCharset =
	std::string("ABCDEFGHIJKLMNOPQRSTUVWXYZ") +
	std::string("abcdefghijklmnopqrstuvwxyz") +
	std::string("1234567890") +
	std::string(".,:;*/+-_#'\"!$%&()[]{}=?<>|~^ ");


FontMap::FontMap() : m_bOk(0)
{
	m_iCharsPerLine = int(std::ceil(std::sqrt(double(m_strCharset.length()))));
	m_iLines = m_iCharsPerLine;

	if(FT_Init_FreeType(&m_ftLib) != 0)
	{
		log_err("Cannot init freetype.");
		return;
	}
}

FontMap::FontMap(const char* pcFont, int iSize) : FontMap()
{
	if(!LoadFont(pcFont, iSize))
	{
		log_err("Cannot load font \"", pcFont, "\".");
		return;
	}
}

bool FontMap::LoadFont(FT_Face ftFace, int iSize)
{
	m_ftFace = ftFace;
	FT_Set_Pixel_Sizes(m_ftFace, 0, iSize);

	m_iTileH = 0;
	m_iTileW = 0;

	int iMaxTop = 0;

	// find /real/ image sizes
	for(std::string::value_type ch : m_strCharset)
	{
		if(FT_Load_Char(m_ftFace, ch,
			FT_RENDER_MODE_NORMAL|FT_LOAD_RENDER /*needed for glyph sizes*/) != 0)
			continue;

		unsigned int iGlyphW = m_ftFace->glyph->bitmap.width;
		unsigned int iGlyphH = m_ftFace->glyph->bitmap.rows;

		FT_Int iTop = m_ftFace->glyph->bitmap_top;
		FT_Pos iRight = m_ftFace->glyph->advance.x >> 6;

		iMaxTop = std::max(int(iTop), iMaxTop);

		m_iTileW = std::max(int(iRight), m_iTileW);
		m_iTileW = std::max(int(iGlyphW), m_iTileW);

		m_iTileH = std::max(int(iGlyphH), m_iTileH);
	}

	m_iPadH = int(nextpow<double>(2, m_iTileH))-m_iTileH;
	m_iPadW = int(nextpow<double>(2, m_iTileW))-m_iTileW;

	m_iLargeW = (m_iTileW+m_iPadW) * m_iCharsPerLine;
	m_iLargeH = (m_iTileH+m_iPadH) * m_iLines;

	// find next larger power of 2 (later used for glTexImage2D)
	m_iLargeW = int(nextpow<double>(2, m_iLargeW));
	m_iLargeH = int(nextpow<double>(2, m_iLargeH));

	m_pcLarge = new unsigned char[m_iLargeW*m_iLargeH];
	memset(m_pcLarge, 0, m_iLargeW*m_iLargeH);

	int iPosY=0, iPosX=0;

	for(std::string::value_type ch : m_strCharset)
	{
		if(FT_Load_Char(m_ftFace, ch, FT_RENDER_MODE_NORMAL|FT_LOAD_RENDER) != 0)
		{
			log_err("Cannot render \"", ch, "\".");
			continue;
		}

		int iGlyphW = m_ftFace->glyph->bitmap.width;
		int iGlyphH = m_ftFace->glyph->bitmap.rows;

		FT_Int iLeft = m_ftFace->glyph->bitmap_left;
		FT_Int iTop = m_ftFace->glyph->bitmap_top;

		const unsigned char* pcBmp = m_ftFace->glyph->bitmap.buffer;
		draw_tile(m_pcLarge, m_iLargeW, m_iLargeH, m_iTileW+m_iPadW, m_iTileH+m_iPadH,
			iPosX, iPosY, pcBmp, iLeft,
			iMaxTop-iTop, iGlyphW, iGlyphH);

		m_mapOffs[ch] = t_offsmap::value_type::second_type(iPosX*(m_iTileW+m_iPadW), iPosY*(m_iTileH+m_iPadH));

		++iPosX;
		if(iPosX >= m_iCharsPerLine)
		{
			++iPosY;
			iPosX = 0;
		}
	}

	m_bOk = 1;
	return true;
}

bool FontMap::LoadFont(const char* pcFont, int iSize)
{
	UnloadFont();

	if(FT_New_Face(m_ftLib, pcFont, 0, &m_ftFace) != 0)
		return false;

	return LoadFont(m_ftFace, iSize);
}

void FontMap::dump(std::ostream& ostr) const
{
	for(unsigned int iY=0; iY<unsigned(m_iLargeH); ++iY)
	{
		for(unsigned int iX=0; iX<unsigned(m_iLargeW); ++iX)
		{
			short sChar = short(m_pcLarge[iY*m_iLargeW + iX]);

			if(sChar)
				ostr << std::hex << std::setw(2) << sChar << " ";
			else
				ostr << "   ";
		}
		ostr << "\n";
	}
}

void FontMap::dump(std::ostream& ostr, const std::pair<int,int>& pair) const
{
	if(pair.first<0 || pair.second<0)
		return;

	for(unsigned int iY=0; iY<unsigned(m_iTileH); ++iY)
	{
		for(unsigned int iX=0; iX<unsigned(m_iTileW); ++iX)
		{
			short sChar = short(m_pcLarge[(iY+pair.second)*m_iLargeW + iX+pair.first]);

			if(sChar)
				ostr << std::hex << std::setw(2) << sChar << " ";
			else
				ostr << "   ";
		}
		ostr << "\n";
	}
}

std::pair<int, int> FontMap::GetOffset(char ch) const
{
	t_offsmap::const_iterator iter = m_mapOffs.find(ch);
	if(iter == m_mapOffs.end())
		return std::pair<int, int>(-1,-1);

	return iter->second;
}


void FontMap::UnloadFont()
{
	m_bOk = 0;

	m_mapOffs.clear();
	if(m_pcLarge) { delete[] m_pcLarge; m_pcLarge = nullptr; }
	if(m_ftFace) { FT_Done_Face(m_ftFace); m_ftFace=0; }
}


FontMap::~FontMap()
{
	UnloadFont();
	if(m_ftLib) { FT_Done_FreeType(m_ftLib); m_ftLib=0; }
}


// ----------------------------------------------------------------------------


template<class T>
GlFontMap<T>::GlFontMap(const char* pcFont, int iSize)
	: FontMap(pcFont, iSize), m_bOk(0)
{
	if(!FontMap::m_bOk)
	{
		log_err("Font map failed to create.");
		return;
	}

	if(!CreateFontTexture())
	{
		log_err("Cannot create font texture.");
		return;
	}

	m_bOk = 1;
}

template<class T>
GlFontMap<T>::GlFontMap(FT_Face ftFace, int iSize)
	: FontMap(), m_bOk(0)
{
	if(!ftFace)
	{
		log_err("Invalid font face.");
		return;
	}

	if(!LoadFont(ftFace, iSize))
	{
		log_err("Font map failed to create.");
		return;
	}

	if(!CreateFontTexture())
	{
		log_err("Cannot create font texture.");
		return;
	}

	m_bOk = 1;
}

template<class T>
bool GlFontMap<T>::CreateFontTexture()
{
	glGetError(); // clear previous errors

	glEnable(GL_TEXTURE_2D);

	glGenTextures(1, &m_tex);
	glBindTexture(GL_TEXTURE_2D, m_tex);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	glTexImage2D(GL_TEXTURE_2D, 0, /*GL_RGBA*/GL_ALPHA, m_iLargeW, m_iLargeH, 0, GL_ALPHA, GL_UNSIGNED_BYTE, m_pcLarge);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR/*_MIPMAP_LINEAR*/ /*GL_NEAREST_MIPMAP_NEAREST*/);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR /*GL_NEAREST*/);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE /*GL_REPEAT*/);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE /*GL_REPEAT*/);

	GLenum glerr = glGetError();
	if(glerr != GL_NO_ERROR)
	{
		log_err("GL error: ", glerr);
		return false;
	}

	return true;
}

template<class T>
void GlFontMap<T>::BindTexture()
{
	if(!m_bOk) return;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, m_tex);

	glEnable(GL_BLEND);
	//glDisable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glDisable(GL_LIGHTING);
	//glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);

	//if(glGetError() != GL_NO_ERROR) log_err("Cannot bind texture.");
}

template<class T>
void GlFontMap<T>::DrawText(T _dX, T _dY, const std::string& str, bool bCenter)
{
	if(!m_bOk) return;

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	T dX = _dX;
	T dY = _dY;

	T dScX = (m_iTileH / 12.) * m_dScale;
	T dScY = T(m_iTileH)/T(m_iTileW) * dScX;

	T dXInc = 1.75*dScX;

	if(bCenter)
	{
		dX -= dXInc * str.length()/2.;
		_dX = dX;
	}

	for(char ch : str)
	{
		if(ch == '\n')
		{
			dX = _dX;
			dY -= 2.*dScY;
			continue;
		}
		else if(ch == '\t')
		{
			ch = ' ';
		}

		std::pair<int, int> offs = GetOffset(ch);
		if(offs.first<0 ||offs.second<0)
		{
			log_err("Character not in map.");
			return;
		}

		T dSizeFullTileW = T(m_iTileW+m_iPadW)/T(m_iLargeW);
		T dSizeFullTileH = T(m_iTileH+m_iPadH)/T(m_iLargeH);

		T dSizeTileW = T(m_iTileW)/T(m_iLargeW);
		T dSizeTileH = T(m_iTileH)/T(m_iLargeH);

		T dPosTileX = T(offs.first/(m_iTileW+m_iPadW))*dSizeFullTileW;
		T dPosTileY = T(offs.second/(m_iTileH+m_iPadH))*dSizeFullTileH;

		glBegin(GL_QUADS);
			gl_traits<T>::SetTextureCoord(dPosTileX, dPosTileY+dSizeTileH);
			gl_traits<T>::SetVertex(-dScX+dX, -dScY+dY);

			gl_traits<T>::SetTextureCoord(dPosTileX+dSizeTileW, dPosTileY+dSizeTileH);
			gl_traits<T>::SetVertex( dScX+dX, -dScY+dY);

			gl_traits<T>::SetTextureCoord(dPosTileX+dSizeTileW, dPosTileY);
			gl_traits<T>::SetVertex( dScX+dX,  dScY+dY);

			gl_traits<T>::SetTextureCoord(dPosTileX, dPosTileY);
			gl_traits<T>::SetVertex(-dScX+dX,  dScY+dY);
		glEnd();

		dX += dXInc;
	}

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}


template<class T>
void GlFontMap<T>::DrawText(T dX, T dY, T dZ, const std::string& str, bool bCenter)
{
	if(!m_bOk) return;

	T dXProj, dYProj;
	gl_proj_pt<t_mat4_gen<T>, t_vec4_gen<T>, T>(dX, dY, dZ, dXProj, dYProj);
	DrawText(dXProj, dYProj, str, bCenter);
}


template<class T>
GlFontMap<T>::~GlFontMap()
{
	glDeleteTextures(GL_TEXTURE_2D, &m_tex);
	m_tex = 0;
}

}
#endif
