source("analysis/map_prep.R")

myv_df %>%
  ggplot(aes(long, lat, group = piece))+
  geom_path()+
  geom_path(data = tibble(long = c(410800, 406900), 
                           lat = c(7277900, 7273500), 
                           piece = 10),
             color = "red", size = 1)+
  coord_equal()

d = tibble(long = c(409600, 406900), 
           lat = c(7278100, 7273500), 
           piece = 10)

myv_df %>%
  filter(piece == 1) %>%
  ggplot(aes(long, lat, group = piece))+
  geom_path()+
  coord_equal()

{myv_df %>%
  filter(piece == 1) %>%
  {pracma::polyarea(x = -.$long, y = .$lat)}}/1000000

test <- myv_df %>%
  filter(piece == 1)


myv_df %>%
  filter(!(long > 409500 & lat > 7278200),
         !(long > 409000 & lat > 7279000)) %>%
  filter(piece == 1) %>%
  ggplot(aes(long, lat, group = piece))+
  geom_path()+
  geom_path(data = tibble(long = c(410800, 406900), 
                          lat = c(7277900, 7273500), 
                          piece = 10),
            color = "red", size = 1)+
  coord_equal()

m = (d$long[2] - d$long[1])/(d$lat[2] - d$lat[1])
b = d$long[2] - m * d$lat[2]

myv_df %>%
  filter(!(long > 409500 & lat > 7278200),
         !(long > 409000 & lat > 7279000),
         long < m * lat  + b) %>%
  ggplot(aes(long, lat, group = piece))+
  geom_path()+
  geom_path(data = tibble(long = c(410800, 406900), 
                          lat = c(7277900, 7273500), 
                          piece = 10),
            color = "red", size = 1)+
  coord_equal()

x

sb <- {myv_df %>%
    filter(!(long > 409500 & lat > 7278200),
           !(long > 409000 & lat > 7279000))  %>%
    filter(piece == 1) %>%
    {pracma::polyarea(x = -.$long, y = .$lat)}}/1000000

wsb <- {myv_df %>%
    filter(!(long > 409500 & lat > 7278200),
           !(long > 409000 & lat > 7279000),
           long < m * lat  + b)  %>%
    filter(piece == 1) %>%
    {pracma::polyarea(x = -.$long, y = .$lat)}}/1000000


wsb/sb * 28.2

17.50917/8.5
