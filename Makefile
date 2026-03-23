# Compilateur
CXX = g++

# Options de compilation
CXXFLAGS = -Wall -Wextra -std=c++11

# Fichiers sources
SRCS = main.cpp Reaction.cpp Entite.cpp

# Fichiers objets
OBJS = $(SRCS:.cpp=.o)

# Exécutable final
TARGET = myFirstGillepsie

# Règle par défaut
all: $(TARGET)
	rm -f $(OBJS)

# Règle pour créer l'exécutable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Règle pour compiler les fichiers .cpp en .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Règle pour nettoyer les fichiers objets et l'exécutable
clean:
	rm -f $(OBJS) $(TARGET)

# Règle pour tout nettoyer et recompiler
rebuild: clean all
