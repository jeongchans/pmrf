#include <gtest/gtest.h>

#include "dpmatrix.h"

TEST(DPMatrixTest, test_dp_matrix) {
    size_t rows = 3;
    size_t cols = 4;
    DPMatrix mat(rows, cols);
    ASSERT_EQ(rows + 1, mat.match.size());
    ASSERT_EQ(rows + 1, mat.insdel.size());
    ASSERT_EQ(rows + 1, mat.delins.size());
    for (size_t i = 0; i <= rows; ++i) {
        EXPECT_EQ(cols + 1, mat.match[i].size());
        EXPECT_EQ(cols + 1, mat.insdel[i].size());
        EXPECT_EQ(cols + 1, mat.delins[i].size());
    }
}

TEST(DPMatrixVisitorTest, test_dp_matrix_visitor) {
    class MyDPMatVisitor : public DPMatrixVisitor {
      public:
        virtual void visit_match(MatchMatrixElement* element, const size_t&, const size_t&) {
            element->score = 0;
            element->from = {NULL, -1, -1};
        }
        virtual void visit_insdel(InsDelMatrixElement* element, const size_t&, const size_t&) {
            element->score = 1;
        }
        virtual void visit_delins(DelInsMatrixElement* element, const size_t&, const size_t&) {
            element->score = 2;
        }
    };

    MyDPMatVisitor visitor;
    size_t rows = 3;
    size_t cols = 4;
    DPMatrix mat(rows, cols);
    mat.traverse(visitor);
    for (size_t i = 0; i <= rows; ++i) {
        for (size_t j = 0; j <= cols; ++j) {
            EXPECT_EQ(0, mat.match[i][j].score);
            EXPECT_EQ(1, mat.insdel[i][j].score);
            EXPECT_EQ(2, mat.delins[i][j].score);
            DPElementPosition prev = mat.match[i][j].from;
            EXPECT_EQ(NULL, prev.element);
            EXPECT_EQ(-1,   prev.row);
            EXPECT_EQ(-1,   prev.col);
        }
    }
}
