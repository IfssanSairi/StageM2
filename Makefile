# Compilateur
CXX = g++

# Options de compilation
CXXFLAGS = -Wall -Wextra -std=c++14

# include paths
INCLUDE = -I/opt/homebrew/Cellar/boost/1.90.0_1/include
#INCLUDE = -I/opt/homebrew/include

# Fichiers sources
SRCS = main.cpp Reaction.cpp Entite.cpp

# Fichiers objets
OBJS = $(SRCS:.cpp=.o)

# ExÃ©cutable final
TARGET = myFirstGillepsie

# RÃ¨gle par dÃ©faut
all: $(TARGET)
	rm -f $(OBJS)

# RÃ¨gle pour crÃ©er l'exÃ©cutable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^
	
# RÃ¨gle pour compiler les fichiers .cpp en .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

# RÃ¨gle pour nettoyer les fichiers objets et l'exÃ©cutable
clean:
	rm -f $(OBJS) $(TARGET)

# RÃ¨gle pour tout nettoyer et recompiler
rebuild: clean all
