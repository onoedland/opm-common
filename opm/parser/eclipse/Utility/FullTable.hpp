/*
  Copyright (C) 2013 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef OPM_PARSER_FULL_TABLE_HPP
#define	OPM_PARSER_FULL_TABLE_HPP

#include "SimpleTable.hpp"

#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cassert>

namespace Opm {
    template <class OuterTable = Opm::SimpleMultiRecordTable, class InnerTable = Opm::SimpleTable>
    class FullTable
    {
        typedef FullTable<OuterTable, InnerTable> Self;
        typedef std::shared_ptr<const OuterTable> OuterTableConstPtr;
        typedef std::shared_ptr<const InnerTable> InnerTableConstPtr;

    protected:
        // protected default constructor for the derived classes
        FullTable() {}

        // protected constructor for the case that the derived classes
        // use specialized classes for the outer and inner tables
        FullTable(Opm::DeckKeywordConstPtr keyword)
        {
            m_outerTable.reset(new OuterTable(keyword));

            for (int rowIdx = 0; rowIdx < static_cast<int>(keyword->size()); ++rowIdx) {
                InnerTableConstPtr curRow(new InnerTable(keyword, /*recordIdx=*/rowIdx));
                m_innerTables.push_back(curRow);
            }
        }

    public:
        typedef std::shared_ptr<Self> Pointer;
        typedef std::shared_ptr<const Self> ConstPointer;

        /*!
         * \brief Read full tables from keywords like PVTO
         *
         * The data for these keywords can be considered a 2D table:
         * The outer one is a multi-record table for a given state,
         * the inner one is a normal table which extends this
         * state. For the PVTO keyword, the outer table represents the
         * gas dissolution factor, pressure, volume factor and
         * viscosity at the oil's saturation point, the inner table is
         * the pressure, volume factor and viscosity of untersaturated
         * oil with the same gas dissolution factor.
         */
        FullTable(Opm::DeckKeywordConstPtr keyword,
                  const std::vector<std::string> &outerColumnNames,
                  const std::vector<std::string> &innerColumnNames)
        {
            m_outerTable.reset(new SimpleMultiRecordTable(keyword, outerColumnNames));

            for (int rowIdx = 0; rowIdx < keyword->size(); ++rowIdx) {
                Opm::SimpleTableConstPtr curRow(
                    new SimpleTable(keyword,
                                    innerColumnNames,
                                    /*recordIdx=*/rowIdx,
                                    /*firstColumnOffset=*/1));
                m_innerTables.push_back(curRow);
            }
        }


        std::shared_ptr<const OuterTable> getOuterTable() const
        { return m_outerTable; }

        std::shared_ptr<const InnerTable> getInnerTable(int rowIdx) const
        {
            assert(0 <= rowIdx && rowIdx < static_cast<int>(m_innerTables.size()));
            return m_innerTables[rowIdx];
        }

    protected:
        std::shared_ptr<const OuterTable> m_outerTable;
        std::vector<std::shared_ptr<const InnerTable> > m_innerTables;

    };

    typedef FullTable<Opm::SimpleMultiRecordTable, Opm::SimpleTable>::Pointer FullTablePtr;
    typedef FullTable<Opm::SimpleMultiRecordTable, Opm::SimpleTable>::ConstPointer FullTableConstPtr;
}

#endif	// OPM_PARSER_FULL_TABLE_HPP

