# Test Set Up
source("../../files.R")
test_source_dir <- 'testcache__'
test_source_subdir <- 'testcache__/subdir'
test_dest_dir <- 'test_data__'
dir.create(test_source_dir)
dir.create(test_source_subdir)
source_txt_top <- sprintf("test%s.txt",1:4)
source_R_top <- sprintf("test%s.R",1:4)
source_txt_sip <- 
new_files <- c(
    file.path(test_source_dir, c(sprintf("test%s.txt",1:4), sprintf("test%s.R", 1:4))),
    file.path(test_source_subdir, "test5.txt"),
)

file.create(new_files, recursive=T)
testthat::test_that("Test files created successfully", {
  testthat::expect_true(all(file.exists(new_files)))
})

# Test
testthat::test_that("Pattern-matched files are moved non-recursively", {
    move_files(from = test_source_dir, to = test_dest_dir, pattern = "[.txt]$", recursive = FALSE)
    expect_out <- file.path(test_dest_dir,sprintf("test%s.txt",1:4))
    testthat::expect_true(all(file.exists(expect_out)))
    testthat::expect_false(file.exists('testdata/txtfiles/test5.txt'))
})

testthat::test_that("Pattern-matched files are moved recursively", {
    move_files(from = test_source_dir, to = "testdata/txtfiles", pattern = "[.txt]$", recursive = TRUE)
    testthat::expect_true(file.exists('testdata/txtfiles/test5.txt'))

})
                    
testthat::test_that("Files do not overwrite", {
    new_in <- file.path(test_source_dir, c(sprintf("test%s.txt",1:4)))
    new_out <- file.path("testdata/txtfiles", c(sprintf("test%s.txt",1:4)))
    file.create(new_in)
    testthat::expect_true(all(file.exists(new_in)), label = "make sure the new test files were generated")
    testthat::expect_true(all(file.exists(new_out)), label = "make sure that expected copy location already contain the same files")
    testthat::expect_message(move_files(from = test_source_dir, to = "testdata/txtfiles", pattern = "[.txt]$", overwrite = FALSE), 
                             regexp='Destination file testdata/txtfiles/test1.txt already exists. File not copied. Resolve copies or set overwrite = TRUE.')
    # attempt to moves the files
    move_files(from = test_source_dir, to = "testdata/txtfiles", pattern = "[.txt]$", overwrite = FALSE)
    testthat::expect_true(all(file.exists(new_in)), label = "new test files still present in source directory")
    testthat::expect_true(all(file.exists(new_out)), label = "'to' file location also contains files of the same name")
})

testthat::test_that("Files do overwrite", {
    new_in <- file.path(test_source_dir, c(sprintf("test%s.txt",1:4)))
    new_out <- file.path("testdata/txtfiles", c(sprintf("test%s.txt",1:4)))
    testthat::expect_true(all(file.exists(new_in)), label = "make sure the new test files were generated")
    testthat::expect_true(all(file.exists(new_out)), label = "make sure that expected copy location already contain the same files")
    # attempt to moves the files
    move_files(from = test_source_dir, to = "testdata/txtfiles", pattern = "[.txt]$", overwrite = TRUE)
    testthat::expect_false(all(file.exists(new_in)), label = "files no longer present in source directory")
    testthat::expect_true(all(file.exists(new_out)), label = "'to' file location contains the files")
})

# test cleanup
system(sprintf('rm -r %s'), test_source_dir)
system(sprintf('rm -r %s'), test_dest_dir)