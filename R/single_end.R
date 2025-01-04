#' Process single-end AB1 files for sequence cleaning.
#'
#' This function takes a directory containing AB1 files, processes them, cleans
#' the sequences based on quality scores, and provides a summary table and a
#' concatenated FASTA output, as well as a list of individual FASTA sequences.
#'
#' @param file_path A character string specifying the path to the directory
#'   containing the AB1 files.
#'
#' @return A list containing:
#'   \itemize{
#'     \item{\code{result_table}}{A data frame summarizing the processing results
#'       for each file, including raw and cleaned sequence lengths, average quality
#'       scores, and deleted base counts.}
#'     \item{\code{fasta_all}}{A character string containing all cleaned sequences
#'       concatenated in FASTA format.}
#'     \item{\code{seq_clean_single_list}}{A list where each element is a cleaned sequence in FASTA format for each corresponding file.}
#'   }
#'
#' @import sangerseqR
#' @import RcppRoll
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' file_path <- "path/to/your/ab1/files"
#' results <- single_end(file_path)
#'
#' # Access results:
#' print(results$result_table)
#' cat(results$fasta_all)
#' print(results$seq_clean_single_list)
#' }
single_end <- function(file_path) {

  if (!dir.exists(file_path)) {
    stop("Error: The specified file path does not exist: ", file_path)
  }

  # 定义一些变量和文件放置的容器
  fasta_all <- character()
  ab1_name <- character()
  result_clean <- data.frame(file_name = character(),
                             length_raw = numeric(),
                             length_clean = numeric(),
                             average_quality_raw = numeric(),
                             average_quality_clean = numeric(),
                             deleted_count = numeric(),
                             deleted_ratio = numeric(),
                             deleted_count1 = numeric(),
                             deleted_count2 = numeric(),
                             deleted_count3 = numeric())
  seq_clean_single_list <- list()

  # 每次运行前把文件夹清空一下
  tempDir <- tempdir()

  seq_clean_single <- file.path(tempDir, "seq_clean_single")

  if (dir.exists(seq_clean_single)) {
    file.remove(list.files(seq_clean_single, full.names = TRUE))
  } else {
    dir.create(seq_clean_single)
  }

  # 取所有文件名
  file_list <- list.files(path = file_path,
                          full.names = TRUE,
                          recursive = FALSE,
                          include.dirs = FALSE
  )

  # 取ab1
  ab1_files <- file_list[endsWith(file_list, ".ab1")]

  # 循环开始
  if (length(ab1_files) > 0){
    for (ab1_file in ab1_files) {

      print(paste("Trying file:", ab1_file)) # 调试

      ab1_name <- tools::file_path_sans_ext(basename(ab1_file)) # 注意这里只取名字没有后缀

      # 读取ab1文件中的碱基和质量信息
      abif_data <- sangerseqR::read.abif(ab1_file)
      quality_values <- abif_data@data$PCON.1
      base_calls <- abif_data@data$PBAS.1
      base_calls_split <- strsplit(base_calls, "")[[1]]
      if (length(quality_values) == length(base_calls_split)) {
        paired_data <- data.frame(base_call = base_calls_split, quality_value = quality_values)
      } else {
        quality_values <- abif_data@data$PCON.2
        base_calls <- abif_data@data$PBAS.2
        base_calls_split <- strsplit(base_calls, "")[[1]]
        paired_data <- data.frame(base_call = base_calls_split, quality_value = quality_values)
      }

      # 进行处理和清洗
      # clean1 用滑窗算法进行降噪和过滤
      window_size <- 10
      quality_averages <- RcppRoll::roll_mean(quality_values, n = window_size, align = "center", fill = 0)

      filter_size <- 45
      windows <- data.frame(num = seq_along(quality_averages), quality_average = quality_averages - filter_size)
      windows

      #将间隔小于10（窗口中心到中心之间距离）的为连续作准备，赋值为0
      transition_indices <- which(diff(sign(windows$quality_average)) != 0) + 1 # 找到断点所在的位置
      transition_indices

      if (length(transition_indices) > 2) {
        for (i in seq(2, length(transition_indices) - 1, by = 2)) { # 找到由正到负的断点
          positive_to_negative_position <- windows$num[transition_indices[i]]
          negative_to_positive_position <- windows$num[transition_indices[i+1]]
          position_diff <- negative_to_positive_position - positive_to_negative_position # 计算gap的距离
          if (position_diff <= 10) { # 用0值填平规定大小的gap，使其被纳入windows filter
            windows$quality_average[positive_to_negative_position:negative_to_positive_position] <- 0
          }
        }
      }

      windows_filter <- subset(windows, quality_average >= 0)
      windows_filter

      if (nrow(windows_filter) == 0) {
        deleted_count1 = nrow(windows)
        deleted_count2 = 0
        deleted_count3 = 0
        paired_data_clean3 = data.frame(base_calls = character(), quality_values = numeric())
      } else {
        deleted_count1 <- nrow(windows) - (nrow(windows_filter) + window_size) # 要把window_size 加上
        deleted_count1
        windows_filter$num

        # clean2 找到最长的连续序列
        x <- windows_filter$num
        current_num <- vector("integer")
        continuous_num <- vector("integer")

        for (i in 1:length(x)) {
          if (i > 1 && x[i] == x[i - 1] + 1) {
            current_num <- c(current_num, x[i])
          } else {
            if (length(current_num) > length(continuous_num)) {
              continuous_num <- current_num
            }
            current_num <- c(x[i])
          }
        }
        if (length(current_num) > length(continuous_num)) {
          continuous_num <- current_num
        }

        deleted_count2 <- length(windows_filter$num) - length(continuous_num)
        paired_data_clean2 <- paired_data[(min(continuous_num) - (window_size/2)):(max(continuous_num) + (window_size/2)), ]

        # clean 3 过滤质量指数小于30的碱基
        threshold <- 30
        paired_data_clean3 <- paired_data_clean2[paired_data_clean2$quality_value >= threshold, ]
        deleted_count3 <- nrow(paired_data_clean2) - nrow(paired_data_clean3)
      }

      # 生成结果表格和对应的fasta文件
      result_row <- data.frame(file_name = ab1_name,
                               length_raw = nrow(paired_data),
                               length_clean = nrow(paired_data_clean3),
                               average_quality_raw = round(mean(paired_data$quality_value), 1),
                               average_quality_clean = round(mean(paired_data_clean3$quality_value), 1),
                               deleted_count = nrow(paired_data) - nrow(paired_data_clean3),
                               deleted_ratio = round((nrow(paired_data) - nrow(paired_data_clean3))/nrow(paired_data), 2),
                               deleted_count1 = deleted_count1,
                               deleted_count2 = deleted_count2,
                               deleted_count3 = deleted_count3)
      result_clean <- rbind(result_clean, result_row)

      fasta_single <- c(">", ab1_name, "\n", paste0(paired_data_clean3$base_call, collapse = ""))

      fasta_all <- c(fasta_all, fasta_single, "\n")

      seq_clean_single_list[[ab1_name]] <- fasta_single

    }

    # 创建输出
    fasta_all <- paste(fasta_all, collapse = " ")
    fasta_all <- gsub(" ", "", fasta_all)

    # 返回结果
    invisible(list(
      result_table = result_clean,
      fasta_all = fasta_all,
      seq_clean_single_list = seq_clean_single_list
    ))

    print("All files were successfully processed!")

  } else {
    message("No .ab1 files found in the specified folder: ", file_path)
  }
}
