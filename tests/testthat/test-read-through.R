test_that(".autonomousDesign() checks its arguments", {
    expect_identical(
        .autonomousDesign(NULL, NULL, 4L),
        matrix(1, 4, 1, dimnames = list(NULL, "(Intercept)"))
    )
    expect_identical(
        ncol(.autonomousDesign(NULL, c("a", "a", "b", "b", "b"), 5L)),
        2L
    )
    expect_error(
        .autonomousDesign(matrix(1, 4, 1), c("a", "b", "a", "b"), 4L),
        "Only one of design and group"
    )
    expect_error(
        .autonomousDesign(NULL, c("a", "b"), 4L),
        "one value for each sample"
    )
    expect_error(
        .autonomousDesign(matrix(1, 3, 1), NULL, 4L),
        "one row for each sample"
    )
    expect_error(
        .autonomousDesign(cbind(1, rep(1, 4)), NULL, 4L),
        "not of full rank"
    )
    expect_error(
        .autonomousDesign(NULL, c("a", "a", "b"), 3L),
        "At least 2 more samples"
    )
})

test_that(".unlistRanges() names ranges by their list element", {
    gr <- GenomicRanges::GRanges(c("chr1:1-10", "chr1:20-30", "chr2:5-8"))
    names(gr) <- c("a", "a", "b")
    grl <- S4Vectors::splitAsList(gr, names(gr))
    expect_identical(names(.unlistRanges(grl)), c("a", "a", "b"))
})
